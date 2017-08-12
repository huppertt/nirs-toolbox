function a = subsasgnDot(a,s,b)
% Assignment to or into a variable.  Any sort of subscripting
% may follow that, and row labels are inherited from the dataset.

%   Copyright 2012-2013 The MathWorks, Inc.


% Translate variable (column) name into an index.
varName = s(1).subs;
if isnumeric(varName) && isscalar(varName)
    % Allow a.(i) where i is an integer
    [varIndex,newNames] = getvarindices(a,varName,true);
    isNewVar = (varIndex > a.nvars);
    if isNewVar
        varName = newNames{1};
    else
        varName = a.varnames{varIndex};
    end
elseif ischar(varName) && size(varName,1) == 1
    varIndex = find(strcmp(varName,a.varnames));
    isNewVar = isempty(varIndex);
    if isNewVar
        % Handle assignment to a property under the 'Properties' (virtual)
        % property, but disallow assignment to properties directly, or to
        % the 'Properties' property.
        if checkreservednames(varName)
            if strcmp(varName,'Properties')
                if ~isscalar(s)
                    try %#ok<ALIGN>
                        a = setproperty(a,s(2:end),b);
                    catch ME, throw(ME); end
                    return
                else
                    error(message('stats:dataset:subsasgn:AssignmentToPropertiesStruct'));
                end
            else % a.ObsNames, a.VarNames
                error(message('stats:dataset:subsasgn:InvalidPropertyAssignment', varName, varName));
            end
        end
        [~, modified] = matlab.lang.makeValidName(varName);
        if any(modified) % error if invalid
            error(message('stats:dataset:InvalidVariableName',varName));
        end
        
        % If this is a new variable, it will go at the end.
        varIndex = a.nvars + 1;
    end
else
    error(message('stats:dataset:subsasgn:IllegalVarSubscript'));
end

% Handle empty assignment intended as deletion of an entire variable or of
% columns/pages/etc. of a variable.  Deletion of rows in a (single)
% variable is caught here and not allowed.  Other empty assignment
% syntaxes may be assignment to cells or may be deletion of things deeper
% in a non-atomic variable, neither is handled here.
if builtin('_isEmptySqrBrktLiteral',b) && (isscalar(s) || isequal(s(2).type,'()'))
    % s(2).type=='()' guarantees that length(s)==2
    if isNewVar
        error(message('stats:dataset:subsasgn:UnrecognizedVarName', varName));
    end
    
    % Syntax:  a.var = []
    %
    % Delete an entire variable.
    if isscalar(s)
        a.data(varIndex) = [];
        a.varnames(varIndex) = [];
        a.nvars = a.nvars - 1;
        % Var-based or properties need to be shrunk.
        if ~isempty(a.props.VarDescription), a.props.VarDescription(varIndex) = []; end
        if ~isempty(a.props.Units), a.props.Units(varIndex) = []; end
        
    % Syntax:  a.var(:,...) = []
    %          a.var(obsIndices,...) = [] is illegal
    %
    % Delete columns/pages/etc. of a variable, with ':' as the first index
    % in subscript.  This may change the dimensionality of the variable,
    % but won't change the number of rows because we require ':' as the
    % first index.
    else
        if ~iscolon(s(2).subs{1})
            error(message('stats:dataset:subsasgn:SingleVariableEmptyAssignmentObs'));
        end
        
        var_j = a.data{varIndex};
        try %#ok<ALIGN>
            var_j(s(2).subs{:}) = [];
        catch ME, throw(ME); end
        a.data{varIndex} = var_j;
    end
    
else
    % Syntax:  a.var = b
    %
    % Replace an entire variable.  It may be shorter than the dataset; it
    % is filled out with default values.  It may be longer than the
    % dataset; existing vars are filled in with default values.  So this
    % is not equivalent to using a colon as the observation index, which
    % cannot change the length of a variable.
    if isscalar(s)
        if isa(b,'dataset')
            error(message('stats:dataset:subsasgn:DatasetVariable'));
        end
        var_j = b;
        newObsNames = {};
        
    % Syntax:  a.var(obsIndices,...) = b
    %          a.var{obsIndices,...} = b
    %          a.var{obsIndices,...} = [] (this is assignment, not deletion)
    %          a.var.field = b
    %
    % Assign to elements in a variable.  Assignment can also be used to
    % expand the variable along a not-first dimension, but expansion
    % operations are not allowed to change the number of rows.
    %
    % Cell indexing, e.g. a.var{obsIndices,...}, or a reference to a
    % field, e.g. a.var.field, may also be followed by deeper levels of
    % subscripting.
    else % ~isscalar(s)
        if isNewVar && (length(s) > 2)
            % Cannot create a new var implicitly by deeper indexing.
            error(message('stats:dataset:subsasgn:UnrecognizedVarName', varName));
        end
        if isequal(s(2).type,'.') % dot indexing into variable
            % No obs labels, but the variable must exist.
            if isNewVar
                error(message('stats:dataset:subsasgn:UnrecognizedVarName', varName));
            end
            var_j = a.data{varIndex};
        else % () or {} subscripting into variable
            % Initialize a new var, or extract an existing var.
            if isNewVar
                % Start the new var out as an empty of the b's class with
                % the same number of rows as the dataset.
                var_j = b(zeros(a.nobs,0));
            else
                var_j = a.data{varIndex};
            end
            
            % The variable inherits observation labels from the dataset.
            % Translate labels to row numbers if necessary.
            obsIndices = s(2).subs{1};
            if iscolon(obsIndices) || islogical(obsIndices) || isnumeric(obsIndices)
                % leave these alone
                newObsNames = {};
            else
                if (size(var_j,2)>1) && isscalar(s(2).subs)
                    error(message('stats:dataset:subsasgn:InvalidLinearIndexing'));
                end
                [obsIndices,~,~,newObsNames] = getobsindices(a, obsIndices, true);
                s(2).subs{1} = obsIndices;
            end
        end
        
        % Now let the variable's subsasgn handle the subscripting in
        % things like a.name(...) or  a.name{...} or a.name.attribute
        
        % Calling subsasgn directly allows changing the shape of var_j by
        % assignment, which the interpreter would not allow.
        if length(s) == 2
            if isobject(var_j)
                if isobject(b) && ~isa(b,class(var_j))
                    try %#ok<ALIGN>
                        var_j = var_j.subsasgn(s(2),b); % dispatch to var_j's subsasgn
                    catch ME, throw(ME); end
                else
                    try %#ok<ALIGN>
                        var_j = subsasgn(var_j,s(2),b);
                    catch ME, throw(ME); end
                end
            else
                % Call builtin, to get correct dispatching even if b is an object.
                try %#ok<ALIGN>
                    var_j = builtin('subsasgn',var_j,s(2),b);
                catch ME, throw(ME); end
            end
        else % length(s) > 2
            % Trick the third and higher levels of subscripting in things like
            % ds.Var{i}(...) etc. into dispatching to the right place even when
            % ds.Var{i}, or something further down the chain, is itself a dataset.
            try %#ok<ALIGN>
                var_j = statslibSubsasgnRecurser(var_j,s(2:end),b);
            catch ME, rethrow(ME); end % point to the line in statslibSubsasgnRecurser
        end
    end
    
    % If this is a new variable, make it official.
    if isNewVar
        a.varnames = [a.varnames varName];
        a.nvars = varIndex;
        % Var-based properties need to be extended.
        if ~isempty(a.props.VarDescription), a.props.VarDescription = [a.props.VarDescription {''}]; end
        if ~isempty(a.props.Units), a.props.Units = [a.props.Units {''}]; end
    end
    
    % If a var was replaced, or a new var was created, and it is
    % shorter than the dataset, fill it out.  It's never the case that
    % assigning into a var can shorten it.
    varLen = size(var_j,1);
    if varLen < a.nobs
        warning(message('stats:dataset:subsasgn:DefaultValuesAddedVariable', varName));
        var_j = lengthenVar(var_j,a.nobs);
    end
    a.data{varIndex} = var_j;
    
    % If a var was expanded by assignment, or if a var was replaced or
    % created, and it is longer than the dataset, fill out the rest of
    % the dataset, including observation names.
    if varLen > a.nobs
        % Don't warn if a had no variables originally
        if a.nvars > 1
            warning(message('stats:dataset:subsasgn:DefaultValuesAdded'));
        end
        a = fillInDataset(a,varLen,newObsNames);
    end
end
