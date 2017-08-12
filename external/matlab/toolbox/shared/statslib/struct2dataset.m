function d = struct2dataset(s,varargin)
%STRUCT2DATASET Convert structure array to dataset array.
%   D = STRUCT2DATASET(S) converts the structure array S to a dataset array D.
%   Each field of S becomes a variable in D.  When S is a scalar structure
%   with N fields, all of which have M rows, then D is an M-by-N array.  When
%   S is a non-scalar M-by-1 structure array with N fields, then D is M-by-N.
%
%   D = STRUCT2DATASET(S, 'PARAM1', VAL1, 'PARAM2', VAL2, ...) specifies optional
%   parameter name/value pairs that determine how the data in S are converted.
%
%      'ReadObsNames'  When ReadObsNames is a field name in S, STRUCT2DATASET uses
%                      that field to create observation names in D, and also sets
%                      D.Properties.DimNames to {ReadObsNames, 'Variables'}.  By
%                      default, or if ReadObsNames is false, STRUCT2DATASET will
%                      not create observation names in D unless you specify names
%                      using ObsNames.
%      'ObsNames'      A cell array of strings containing observation names for D.
%                      The names need not be valid MATLAB identifiers, but must be
%                      unique.
%      'AsScalar'      A logical value indicating whether or not to treat a scalar
%                      structure S specially, as described above.  Setting this to
%                      false causes STRUCT2DATASET to convert S to a dataset with
%                      SIZE(S,1) observations, even when S is a scalar structure.
%                      Default is true when S is a scalar structure, false otherwise.
%
%   See also DATASET2STRUCT, CELL2DATASET, DATASET.

%   Copyright 2012 The MathWorks, Inc.


if ~isstruct(s)
    error(message('stats:struct2dataset:NotColumn'));
end

pnames = {'ObsNames' 'ReadObsNames' 'AsScalar'};
dflts =  {       []          false         [] };
[obsnames,readObsnames,asScalar,supplied] ...
    = statslib.internal.parseArgs(pnames, dflts, varargin{:});

if supplied.AsScalar
    if asScalar && ~isscalar(s)
        error(message('stats:struct2dataset:NonScalar'));
    end
else
    asScalar = isscalar(s);
end

haveObsnames = false;
if supplied.ObsNames
    if readObsnames
        error(message('stats:struct2dataset:ObsNamesParamConflict'));
    end
    haveObsnames = true;
elseif readObsnames
    if ischar(readObsnames)
        obsnamesField = readObsnames;
    else % true, read from the first field
        obsnamesField = fieldnames(s); obsnamesField = obsnamesField{1};
    end
    if asScalar
        obsnames = s.(obsnamesField);
    else
        obsnames = {s.(obsnamesField)}';
    end
    s = rmfield(s,obsnamesField);
    haveObsnames = true;
end

if asScalar
    if isempty(fieldnames(s))
        % Size the array according to the obs names
        d = dataset('ObsNames',obsnames);
    else
        d = dataset(s);
        if haveObsnames, d.Properties.ObsNames = obsnames; end
    end
else
    if ~iscolumn(s)
        error(message('stats:struct2dataset:NotColumn'));
    end
    
    fnames = fieldnames(s);
    nvars = length(fnames);
    vars = cell(1,nvars);
    for j = 1:nvars
        cj = {s.(fnames{j})}';
        if isempty(cj) % prevent iscellstr from catching these
            % give these the right number of rows, but no columns
            vars{j} = zeros(size(cj,1),0);
        elseif iscellstr(cj)
            % Prevent a cellstr that happens to have all the same length strings,
            % e.g., datestrs, from being converted into a char matrix.
            vars{j} = cj;
        elseif any(cellfun(@(x)size(x,1),cj(:)) ~= 1)
            % If the cells don't all have one row, we won't be able to
            % concatenate them and preserve observations, leave it as is.
            vars{j} = cj;
        else
            % Concatenate cell contents into a homogeneous array (if all cells
            % of cj contain "atomic" values), a cell array (if all cells of cj
            % contain cells), or an object array (if all the cells of cj
            % contain objects).  The result may have multiple columns or pages
            % if the cell contents did, but each row will correspond to a
            % "row" (i.e., element) of S.  If that fails, leave it as a cell.
            try
                vars_j = cell(1,size(cj,2));
                % Concatenate rows first
                for i = 1:size(cj,2), vars_j{i} = cat(1,cj{:,i}); end
                % Now concatenate multiple columns into a matrix
                vars{j} = cat(2,vars_j{:});
            catch ME %#ok<NASGU>
                vars{j} = cj;
            end
        end
    end
    numRows = numel(s);
    if isempty(vars) % creating a dataset with no variables
        % Give the output dataset the same number of rows as the input struct ...
        if haveObsnames % ... using either the supplied observation names
            d = dataset('ObsNames',obsnames);
        else            % ... or by tricking the constructor
            dummyNames = cellstr(num2str((1:numRows)'));
            d = dataset('ObsNames',dummyNames(1:numRows));
            d.Properties.ObsNames = {};
        end
    else
        d = dataset(cell2struct(vars,fnames,2));
        if haveObsnames, d.Properties.ObsNames = obsnames; end
    end
end

if readObsnames
    d.Properties.DimNames{1} = obsnamesField;
end
