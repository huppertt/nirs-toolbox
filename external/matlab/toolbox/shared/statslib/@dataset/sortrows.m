function [b,idx] = sortrows(a,vars,sortMode)
%SORTROWS Sort rows of a dataset array.
%   B = SORTROWS(A) returns a copy of the dataset A, with the observations
%   sorted in ascending order by all of the variables in A.  The observations
%   in B are sorted first by the first variable, next by the second variable,
%   and so on.  Each variable in A must be a valid input to SORT, or, if the
%   variable has multiple columns, to the MATLAB SORTROWS function or to its
%   own SORTROWS method.
%
%   B = SORTROWS(A,VARS) sorts the observations in A by the variables
%   specified by VARS.  VARS is a positive integer, a vector of positive
%   integers, a variable name, a cell array containing one or more variable
%   names, or a logical vector.
%
%   B = SORTROWS(A,'ObsNames') sorts the observations in A by the observation
%   names.
%
%   B = SORTROWS(A,VARS,MODE) sorts A in the direction(s) specified by MODE.
%   When MODE is 'ascend' (the default) or 'descend', SORTROWS sorts A in
%   ascending or descending order, respectively, for all variables specified
%   by VARS.  MODE may also be a cell array containing the strings 'ascend' or
%   'descend' to specify a different direction for each variable specified by
%   VARS.  Specify VARS as [] to sort using all variables.
%
%   [B,IDX] = SORTROWS(A, ...) also returns an index vector IDX such that
%   B = A(IDX,:).
%
%   See also DATASET/UNIQUE.

%   Copyright 2006-2012 The MathWorks, Inc.


if nargin < 2 || isempty(vars)
    vars = 1:a.nvars;
elseif strcmpi(vars,'obsnames')
    % When sorting on the observation names, no other variables are allowed,
    % and there'd be no point since the names are unique
    vars = 0; % special flag
else % translate names or logical to indices
    vars = getvarindices(a,vars,false);
end

sortModeStrs = {'ascend','descend'};
if nargin < 3 || isempty(sortMode)
    sortMode = ones(size(vars));
else
    if ischar(sortMode)
        sortMode = cellstr(sortMode);
    elseif ~iscellstr(sortMode)
        error(message('stats:dataset:sortrows:UnrecognizedMode'));
    end
    [tf,sortMode] = ismember(lower(sortMode(:)),sortModeStrs); % 1 or 2
    if ~all(tf)
        error(message('stats:dataset:sortrows:UnrecognizedMode'));
    elseif isscalar(sortMode)
        sortMode = repmat(sortMode,size(vars));
    elseif length(sortMode) ~= length(vars)
        error(message('stats:dataset:sortrows:WrongLengthMode'));
    end
end

% Sort on each index variable, last to first.  Since sort is stable, the
% result is as if they were sorted all together.
if vars == 0 % 'ObsNames'
    if isempty(a.obsnames)
        warning(message('stats:dataset:sortrows:EmptyObsNames'));
    else
        [~,idx] = sort(a.obsnames);
        if sortMode==2, idx = flipud(idx); end % obsnames are unique => stable
    end
else
    idx = (1:a.nobs)';
    a_data = a.data;
    for j = length(vars):-1:1
        var_j = a_data{vars(j)};
        if ~ismatrix(var_j)
            error(message('stats:dataset:sortrows:NDVar',a.varnames{vars(j)}));
        end
        var_j = var_j(idx,:);
        % cell/sort is only for cellstr, use sortrows for cell always.
        if ~iscell(var_j) && isvector(var_j) && (size(var_j,2) == 1)
            try
                [~,ord] = sort(var_j,1,sortModeStrs{sortMode(j)});
            catch ME
                error(message('stats:dataset:sortrows:SortOnVarFailed',a.varnames{vars(j)},class(var_j)));
            end
        else % multi-column, or cell
            % Sort by all columns, all either ascending or descending
            cols = (1:size(var_j,2)) * 2*(1.5-sortMode(j));
            try
                [~,ord] = sortrows(var_j,cols);
            catch ME
                if iscell(var_j) && ~all(cellfun(@isscalar,var_j(:)))
                    error(message('stats:dataset:sortrows:SortOnNonScalarCellValues',a.varnames{vars(j)}));
                else
                    error(message('stats:dataset:sortrows:SortrowsOnVarFailed',a.varnames{vars(j)},class(var_j)));
                end
            end
        end
        idx = idx(ord);
    end
end

s.type = '()'; s.subs = {idx ':'};
b = subsrefParens(a,s);
