function [b,idx] = sortrows(a,vars,sortMode)
%SORTROWS Sort rows of a table.
%   B = SORTROWS(A) returns a copy of the table A, with the rows sorted in
%   ascending order by all of the variables in A.  The rows in B are sorted
%   first by the first variable, next by the second variable, and so on.  Each
%   variable in A must be a valid input to SORT, or, if the variable has
%   multiple columns, to the MATLAB SORTROWS function or to its own SORTROWS
%   method.
%
%   B = SORTROWS(A,VARS) sorts the rows in A by the variables specified by VARS.
%   VARS is a positive integer, a vector of positive integers, a variable name,
%   a cell array containing one or more variable names, or a logical vector.
%
%   VARS may also contain a mix of positive and negative integers.  If an
%   element of VARS is positive, the corresponding variable in A will be sorted
%   in ascending order; if an element of VARS is negative, the corresponding
%   variable in A will be sorted in descending order.  These signs are ignored
%   if you provide the MODE input described below.
%
%   B = SORTROWS(A,'RowNames') sorts the rows in A by the row names.
%
%   B = SORTROWS(A,VARS,MODE) sorts A in the direction(s) specified by MODE.
%   When MODE is 'ascend' (the default) or 'descend', SORTROWS sorts A in
%   ascending or descending order, respectively, for all variables specified
%   by VARS.  MODE may also be a cell array containing the strings 'ascend' or
%   'descend' to specify a different direction for each variable specified by
%   VARS.  Specify VARS as 1:SIZE(A,2) to sort using all variables.
%
%   [B,IDX] = SORTROWS(A, ...) also returns an index vector IDX such that
%   B = A(IDX,:).
%
%   See also UNIQUE.

%   Copyright 2012-2013 The MathWorks, Inc.

sortSigns = [];
if nargin < 2 % do not treat [] as "default behavior"
    vars = 1:a.nvars;
elseif strcmp(vars,'RowNames')
    % When sorting on the row names, no other variables are allowed,
    % and there'd be no point since the names are unique
    vars = 0; % special flag
else % translate names or logical to indices
    if isnumeric(vars)
        if nargin < 3 || isempty(sortMode)
            sortSigns = sign(vars);
        end
        vars = abs(vars);
    end
    vars = getVarIndices(a,vars,false);
end

sortModeStrs = {'ascend','descend'};
if nargin < 3 || isempty(sortMode)
    if isempty(sortSigns)
        sortMode = ones(size(vars));
    else
        sortMode = 1 + (sortSigns == -1); % 1 or 2
    end
else
    if ischar(sortMode)
        sortMode = cellstr(sortMode);
    elseif ~iscellstr(sortMode)
        error(message('MATLAB:table:sortrows:UnrecognizedMode'));
    end
    [tf,sortMode] = ismember(lower(sortMode(:)),sortModeStrs); % 1 or 2
    if ~all(tf)
        error(message('MATLAB:table:sortrows:UnrecognizedMode'));
    elseif isscalar(sortMode)
        sortMode = repmat(sortMode,size(vars));
    elseif length(sortMode) ~= length(vars)
        error(message('MATLAB:table:sortrows:WrongLengthMode'));
    end
end

% Sort on each index variable, last to first.  Since sort is stable, the
% result is as if they were sorted all together.
if vars == 0 % 'RowNames'
    if isempty(a.rownames)
        warning(message('MATLAB:table:sortrows:EmptyRowNames'));
        idx = 1:a.nrows;
    else
        [~,idx] = sort(a.rownames);
        if sortMode==2, idx = flipud(idx); end % rownames are unique => stable
    end
else
    idx = (1:a.nrows)';
    a_data = a.data;
    for j = length(vars):-1:1
        try
            var_j = a_data{vars(j)};
        catch ME
            if any(strcmp(ME.identifier,{'MATLAB:badsubscript' 'MATLAB:cellIndexIsZero'}))
                % Non-numeric var indices have already been converted to numeric, an
                % error here is from a numeric index that was bad to begin with
                error(message('MATLAB:table:sortrows:BadNumericVarIndices'));
            else
                rethrow(ME);
            end
        end
        if ~ismatrix(var_j)
            error(message('MATLAB:table:sortrows:NDVar',a.varnames{vars(j)}));
        end
        var_j = var_j(idx,:);
        % cell/sort is only for cellstr, use sortrows for cell always.
        if ~iscell(var_j) && isvector(var_j) && (size(var_j,2) == 1)
            try
                [~,ord] = sort(var_j,1,sortModeStrs{sortMode(j)});
            catch ME
                m = message('MATLAB:table:sortrows:SortOnVarFailed',a.varnames{vars(j)},class(var_j));
                throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
            end
        else % multi-column, or cell
            % Sort by all columns, all either ascending or descending
            cols = (1:size(var_j,2)) * 2*(1.5-sortMode(j));
            try
                [~,ord] = sortrows(var_j,cols);
            catch ME
                m = message('MATLAB:table:sortrows:SortrowsOnVarFailed',a.varnames{vars(j)},class(var_j));
                throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
            end
        end
        idx = idx(ord);
    end
end

b = subsrefParens(a,{idx ':'});
