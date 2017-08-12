function c = dataset2cell(d,varargin)
%DATASET2CELL Convert dataset array to cell array.
%   C = DATASET2CELL(D) converts the dataset array D to a cell array C.  Each
%   variable of D becomes a column in C.  If D is an M-by-N array, then C is
%   (M+1)-by-N, with D's variable names in the first row.  If D contains
%   observation names, then C is (M+1)-by-(N+1), with those names in the first
%   column.
%
%   C = DATASET2CELL(D, 'PARAM1', VAL1, 'PARAM2', VAL2, ...) specifies optional
%   parameter name/value pairs that determine how the data in D are converted.
%
%      'SplitVars'  A logical value indicating whether variables in D that
%                   have multiple columns should be split apart to become
%                   multiple columns in C, or be kept intact as a single
%                   column in C.  The default is false.
%
%   See also CELL2DATASET, DATASET2STRUCT, DATASET.

%   Copyright 2011-2012 The MathWorks, Inc.


pnames = {'SplitVars'};
dflts =  {     false };
[splitVars] = statslib.internal.parseArgs(pnames, dflts, varargin{:});

[nobs,nvars] = size(d);
haveObsNames = ~isempty(d.obsnames);
if splitVars
    % Each column of each variable in D becomes a column in C, with a suitable
    % suffix for the column heading.  N-D variables are coerced to 2D.  Create
    % a mapping from vars to cell columns.
    varNumCols = datasetfun(@(x)size(x(:,:),2),d);
    c = cell(nobs+1,sum(varNumCols)+haveObsNames);
    var2cols = cell(1,length(varNumCols));
    colNames = cell(1,sum(varNumCols));
    colsj = 0;
    for j = 1:length(var2cols)
        if varNumCols(j) > 0
            colsj = colsj(end) + (1:varNumCols(j));
            var2cols{j} = colsj + haveObsNames;
            colNames(colsj) = makeColNames(d.varnames(j),varNumCols(j));
        end
    end
else
    % Each variable in D becomes a single column in C.
    c = cell(nobs+1,nvars+haveObsNames);
    var2cols = num2cell((1:size(d,2)) + haveObsNames);
    colNames = d.varnames;
end

if haveObsNames
    c{1,1} = 'ObsNames';
    c(2:end,1) = d.obsnames;
    c(1,2:end) = colNames;
else
    c(1,:) = colNames;
end

d_data = d.data;
for j = 1:nvars
    vj = d_data{j};
    if iscell(vj)
        if iscolumn(vj) || splitVars
            % If the cell var is a single column, or if we're splitting up
            % vars, then copy it into D as is.
            c(2:end,var2cols{j}) = vj;
        else
            % If the cell var is multi-column and we're not splitting up vars,
            % then break it apart by rows, but keep each row intact.
            c(2:end,var2cols{j}) = mat2cell(vj,ones(nobs,1));
        end
    else
        % If the variable is not a cell array, split it up into cells, one per
        % row, and perhaps one per column.
        if splitVars
            if isnumeric(vj)
                c(2:end,var2cols{j}) = num2cell(vj(:,:));
            else % non-numeric, but not cell
                c(2:end,var2cols{j}) = mat2cell(vj(:,:),ones(1,nobs),ones(1,varNumCols(j))); %#ok<MMTC>
            end
        else
            c(2:end,var2cols{j}) = mat2cell(vj,ones(nobs,1));
        end
    end
end


function names = makeColNames(baseName,n)
if n == 1
    names = baseName;
else
    names = strcat(baseName,num2str((1:n)','_%-d'));
end
