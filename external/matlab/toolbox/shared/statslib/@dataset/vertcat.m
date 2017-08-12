function a = vertcat(varargin)
%VERTCAT Vertical concatenation for dataset arrays.
%   DS = VERTCAT(DS1, DS2, ...) vertically concatenates the dataset arrays
%   DS1, DS2, ... .  Observation names, when present, must be unique across
%   datasets.  VERTCAT fills in default observation names for the output when
%   some of the inputs have names and some do not.
%
%   Variable names for all dataset arrays must be identical except for order.
%   VERTCAT concatenates by matching variable names.  VERTCAT assigns values
%   for the "per-variable" properties (e.g., Units and VarDescription) in DS
%   from the corresponding property values in DS1.
%
%   See also DATASET/CAT, DATASET/HORZCAT.

%   Copyright 2006-2012 The MathWorks, Inc.


b = varargin{1};
if isequal(b,[]) % accept this as a valid "identity element"
    b = dataset;
elseif ~isa(b,'dataset')
    error(message('stats:dataset:vertcat:InvalidInput'));
end
a = b;
[a_varnames,a_varord] = sort(b.varnames);
for j = 2:nargin
    b = varargin{j};
    if isequal(b,[]) % accept this as a valid "identity element"
        b = dataset;
    elseif ~isa(b,'dataset')
        error(message('stats:dataset:vertcat:InvalidInput'));
    end
    
    % some special cases to mimic built-in behavior
    if a.nvars==0 && a.nobs==0
        a = b;
        [a_varnames,a_varord] = sort(b.varnames);
        continue;
    elseif b.nvars==0 && b.nobs==0
        % do nothing
        continue;
    elseif a.nvars ~= b.nvars
        error(message('stats:dataset:vertcat:SizeMismatch'));
    end
    
    [b_varnames,b_varord] = sort(b.varnames);
    if ~all(strcmp(a_varnames,b_varnames))
        error(message('stats:dataset:vertcat:UnequalVarNames'));
    end
    
    if ~isempty(a.obsnames) && ~isempty(b.obsnames)
        checkduplicatenames(b.obsnames,a.obsnames,'obsnames');
        a.obsnames = vertcat(a.obsnames, b.obsnames);
    elseif ~isempty(b.obsnames) % && isempty(a.obsnames)
        a_obsnames = dfltobsnames(1:a.nobs);
        b_obsnames = matlab.lang.makeUniqueStrings(b.obsnames, a_obsnames, namelengthmax);
        a.obsnames = vertcat(a_obsnames, b_obsnames);
    elseif ~isempty(a.obsnames) % && isempty(b.obsnames)        
        b_obsnames = dfltobsnames(a.nobs+(1:b.nobs));
        b_obsnames = matlab.lang.makeUniqueStrings(b_obsnames, a.obsnames, namelengthmax);
        a.obsnames = vertcat(a.obsnames, b_obsnames);
    end

    b_reord(a_varord) = b_varord; %#ok<AGROW>, full reassignment each time
    a.nobs = a.nobs + b.nobs;
    for i = 1:a.nvars
        % Prevent concatenation of a cell variable with a non-cell variable, which
        % would add only a single cell to the former, containing the latter.
        if iscell(a.data{i}) && ~iscell(b.data{b_reord(i)})
            error(message('stats:dataset:vertcat:VertcatCellAndNonCell', a.varnames{ i }));
        end
        try
            a.data{i} = vertcat(a.data{i}, b.data{b_reord(i)});
        catch ME
            m = message('stats:dataset:vertcat:VertcatMethodFailed',a.varnames{i});
            throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
        end
        % Something is badly wrong with whatever vertcat method has been called.
        if size(a.data{i},1) ~= a.nobs
            error(message('stats:dataset:vertcat:VertcatWrongLength', a.varnames{ i }));
        end
    end
end
