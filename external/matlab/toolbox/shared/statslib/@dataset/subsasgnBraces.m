function a = subsasgnBraces(a,s,b)
% '{}' is assignment to or into the contents of a subset of a dataset array.
% Any sort of subscripting may follow.

%   Copyright 2012-2013 The MathWorks, Inc.


% '{}' is assignment of raw values into a dataset element.  Could be any
% sort of subscript following that.  The shape of the element differs,
% depending on the dimensionality of the var: if the var is nxp, the
% element is 1xp, while if the var is nxpxqx..., the element is pxqx... .
% This is much like time series behavior.  Also, if the var is a column of
% cells, then the element is technically a scalar cell, but it seems
% sensible to do one extra "contents of", and not force callers to say
% a{i,j}{1}.
if numel(s(1).subs) ~= a.ndims
    error(message('stats:dataset:subsasgn:NDSubscript'));
end

% Translate observation (row) names into indices (leaves ':' alone)
[obsIndex,numObsIndices,~,newObsNames] = getobsindices(a, s(1).subs{1}, true);

% Translate variable (column) names into indices (translates ':').  Do not
% allow variable creation with {}-indexing.
varIndex = getvarindices(a, s(1).subs{2}, false);

if numObsIndices > 1 || ~isscalar(varIndex)
    error(message('stats:dataset:subsasgn:MultipleElementAssignment'));
end

% Extract an existing var
var_j = a.data{varIndex};

% Syntax:  a{obsIndex,varIndex} = b
%
% Assignment to an element of a dataset.
if isscalar(s)
    if builtin('_isEmptySqrBrktLiteral',b) && ~iscell(var_j)
        error(message('stats:dataset:subsasgn:InvalidEmptyAssignmentToElement'));
    elseif iscell(var_j)
        if numel(var_j) == size(var_j,1)
            % If the element is a scalar cell, assign into its contents.
            var_j{obsIndex,:} = b;
        else
            error(message('stats:dataset:subsasgn:MultipleCellAssignment'));
        end
    else
        % Set up a subscript expression that will assign to the entire
        % element for the specified observation/variable.  Size checks
        % will be handled by a{i,j}'s subsasgn.
        subs{1} = obsIndex; subs{2:ndims(var_j)} = ':';
        try %#ok<ALIGN>
            var_j(subs{:}) = b;
        catch ME, throw(ME); end
        % *** this error may not even be possible ***
        if size(var_j,1) ~= a.nobs
            error(message('stats:dataset:subsasgn:InvalidVarReshape'));
        end
    end

% Syntax:  a{obsIndex,varIndex}(...) = b
%          a{obsIndex,varIndex}{...} = b
%          a{obsIndex,varIndex}.name = b
%
% Assignment into an element of a dataset.  This operation is allowed
% to change the shape of the variable, as long as the number of rows
% does not change.
else % ~isscalar(s)
    if iscell(var_j)
        if numel(var_j) == size(var_j,1)
            % If the element is a scalar cell, assign into its contents
            s(1).subs = {obsIndex}; % s(1).type is already '{}'
        else
            error(message('stats:dataset:subsasgn:MultipleCellAssignment'));
        end

    else
        % Transfer the observation index from the dataset-level
        % subscript expression to the beginning of the existing
        % element subscript expression, and do the assignment at
        % the element level.
        s(2).subs = [obsIndex s(2).subs];
        s = s(2:end);
    end

    % Let a{i,j}'s subsasgn handle the cascaded subscript expressions.

    % *** subsasgn allows certain operations that the interpreter
    % *** would not, for example, changing the shape of var_j by
    % *** assignment.
    if isscalar(s) % ~iscell(var_j) && length(s_original)==2
        if isobject(var_j)
            if isobject(b) && ~isa(b,class(var_j))
                try %#ok<ALIGN>
                    var_j = var_j.subsasgn(s,b); % dispatch to var_j's subsasgn
                catch ME, throw(ME); end
            else
                try %#ok<ALIGN>
                    var_j = subsasgn(var_j,s,b);
                catch ME, throw(ME); end
            end
        else
            % Call builtin, to get correct dispatching even if b is an object.
            try %#ok<ALIGN>
                var_j = builtin('subsasgn',var_j,s,b);
            catch ME, throw(ME); end
        end
    else % ~iscell(var_j) && length(s_original)>2, or iscell(var_j) && length(s_original)>1
        % *** A hack to get the third and higher levels of subscripting in
        % *** things like ds{i,'Var'}(...) etc. to dispatch to the right place
        % *** when ds{i,'Var'}, or something further down the chain, is itself
        % *** a dataset.
        try %#ok<ALIGN>
            var_j = statslibSubsasgnRecurser(var_j,s,b);
        catch ME, rethrow(ME); end % point to the line in statslibSubsasgnRecurser
    end

    % Do not allow growing a variable with brace assignment into a variable
    if size(var_j,1) ~= a.nobs
        error(message('stats:dataset:subsasgn:InvalidVarReshape'));
    end
end

% If the var is shorter than the dataset, fill it out.  This should never
% happen; assigning into a var cannot shorten the number of rows.
varLen = size(var_j,1);
if varLen < a.nobs
    warning(message('stats:dataset:subsasgn:DefaultValuesAddedVariable', a.varnames{varIndex}));
    var_j = lengthenVar(var_j, a.nobs);

% If a var was lengthened by assignment, fill out the rest of the dataset,
% including observation names.
elseif varLen > a.nobs
    warning(message('stats:dataset:subsasgn:DefaultValuesAdded'));
    a = fillInDataset(a,varLen,newObsNames);
end

a.data{varIndex} = var_j;
