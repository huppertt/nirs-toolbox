function a = subsasgnParens(a,s,b,creating)
% '()' is assignment to a subset of a dataset array.  Only dot subscripting
% may follow.

%   Copyright 2012-2013 The MathWorks, Inc.


if numel(s(1).subs) ~= a.ndims
    error(message('stats:dataset:subsasgn:NDSubscript'));
end

if ~isscalar(s)
    switch s(2).type
    case '()'
        error(message('stats:dataset:subsasgn:InvalidSubscriptExpr'));
    case '{}'
        error(message('stats:dataset:subsasgn:InvalidSubscriptExpr'));
    case '.'
        if creating
            error(message('stats:dataset:subsasgn:InvalidSubscriptExpr'));
        end
        
        % Syntax:  a(obsIndices,varIndices).name = b
        %
        % Assignment into a variable of a subarray.
        
        % Get the subarray, do the dot-variable assignment on that
        % *** need to accept out of range obs or var indices/names here
        % *** as assignment to grow the dataset
        c = subsrefParens(a,s(1));
        b = subsasgnDot(c,s(2:end),b);
        
        % Now let the simple () subscripting code handle assignment of the updated
        % subarray back into the original array.
        s = s(1);
    end
end


% If a new dataset is being created, or if the LHS is 0x0, then interpret
% ':' as the size of the corresponding dim from the RHS, not as nothing.
deleting = builtin('_isEmptySqrBrktLiteral',b);
colonFromRHS = ~deleting && (creating || all(size(a)==0));

% Translate observation (row) names into indices (leave ':' alone)
if colonFromRHS && iscolon(s(1).subs{1})
    obsIndices = 1:b.nobs;
    numObsIndices = b.nobs;
    maxObsIndex = b.nobs;
    newObsNames = {};
else
    [obsIndices,numObsIndices,maxObsIndex,newObsNames] = ...
                               getobsindices(a, s(1).subs{1}, ~deleting);
end

% Translate variable (column) names into indices (translate ':' to 1:nvars)
if colonFromRHS && iscolon(s(1).subs{2})
    varIndices = 1:b.nvars;
    newVarNames = b.varnames;
else
    [varIndices,newVarNames] = getvarindices(a, s(1).subs{2}, ~deleting);
end

% Syntax:  a(obsIndices,:) = []
%          a(:,varIndices) = []
%          a(obsIndices,varIndices) = [] is illegal
%
% Deletion of complete observations or entire variables.
if deleting
    % Delete observations across all variables
    if iscolon(s(1).subs{2})
        if isnumeric(obsIndices)
            obsIndices = unique(obsIndices);
            numObsIndices = numel(obsIndices);
        end
        newNobs = a.nobs - numObsIndices;
        a_data = a.data;
        for j = 1:a.nvars
            var_j = a_data{j};
            if ismatrix(var_j)
                var_j(obsIndices,:) = []; % without using reshape, may not be one
            else
                sizeOut = size(var_j); sizeOut(1) = newNobs;
                var_j(obsIndices,:) = [];
                var_j = reshape(var_j,sizeOut);
            end
            a_data{j} = var_j;
        end
        a.data = a_data;
        if ~isempty(a.obsnames), a.obsnames(obsIndices) = []; end
        a.nobs = newNobs;

        % Delete entire variables
    elseif iscolon(s(1).subs{1})
        varIndices = unique(varIndices); % getvarindices converts all varindex types to numeric
        a.data(varIndices) = [];
        a.varnames(varIndices) = [];
        a.nvars = a.nvars - numel(varIndices);
        % Var-based properties need to be shrunk.
        if ~isempty(a.props.VarDescription), a.props.VarDescription(varIndices) = []; end
        if ~isempty(a.props.Units), a.props.Units(varIndices) = []; end

    else
        error(message('stats:dataset:subsasgn:InvalidEmptyAssignment'));
    end

% Syntax:  a(obsIndices,varIndices) = b
%
% Assignment from a dataset.  This operation is supposed to replace or
% grow at the level of the _dataset_.  So no internal reshaping of
% variables is allowed -- we strictly enforce sizes. In other words, the
% existing dataset has a specific size/shape for each variable, and
% assignment at this level must respect that.
elseif isa(b,'dataset')
    if isscalar(b) % scalar expansion of a single dataset element, which may itself be non-scalar
        b = scalarRepmat(b,numObsIndices,length(varIndices));
    else
        if b.nobs ~= numObsIndices
            error(message('stats:dataset:subsasgn:ObsDimensionMismatch'));
        end
        if b.nvars ~= length(varIndices)
            error(message('stats:dataset:subsasgn:VarDimensionMismatch'));
        end
    end

    existingVarLocs = find(varIndices <= a.nvars);
    a_data = a.data; b_data = b.data;
    for j = existingVarLocs
        var_j = a_data{varIndices(j)};
        % The size of the RHS has to match what it's going into.
        sizeLHS = size(var_j); sizeLHS(1) = numObsIndices;
        if ~isequal(sizeLHS, size(b_data{j}))
            error(message('stats:dataset:subsasgn:DimensionMismatch', a.varnames{varIndices(j)}));
        end
        if iscolon(obsIndices)
            var_j = b_data{j};
        else
            try %#ok<ALIGN>
                var_j(obsIndices,:) = b_data{j}(:,:);
            catch ME, throw(ME); end
        end
        % No need to check for size change, RHS and LHS are identical sizes.
        a_data{varIndices(j)} = var_j;
    end
    a.data = a_data;

    % Add new variables if necessary.  Note that b's varnames do not
    % propagate to a in () assignment, unless a is being created or grown
    % from 0x0.  They do for horzcat, though.
    newVarLocs = find(varIndices > a.nvars);
    if ~isempty(newVarLocs)
        [~, modified] = matlab.lang.makeValidName(newVarNames);
        if any(modified) % error if invalid
            firstModified = newVarNames{find(modified,1)};
            error(message('stats:dataset:InvalidVariableName',firstModified));
        end
        
        a_data = [a_data cell(1,length(newVarNames))];
        a_nobs = a.nobs;
        for j = 1:length(newVarNames)
            var_b = b_data{newVarLocs(j)};
            if iscolon(obsIndices)
                var_j = var_b;
            else
                % Start the new variable out as 0-by-(trailing size of b),
                % then let the assignment add rows.
                var_j = repmat(var_b,[0 ones(1,ndims(var_b)-1)]);
                var_j(obsIndices,:) = var_b(:,:);
            end
            % A new var may need to grow to fit the dataset
            if size(var_j,1) < a_nobs
                warning(message('stats:dataset:subsasgn:DefaultValuesAddedVariable', newVarNames{j}));
                var_j = lengthenVar(var_j, a_nobs);
            end
            a_data{a.nvars+j} = var_j;
        end
        a.data = a_data;
        LHSVars = 1:a.nvars;
        RHSVars = 1:b.nvars;
        a.varnames = [a.varnames newVarNames];
        a.nvars = a.nvars + length(newVarNames);
        % Var-based properties need to be extended.
        a.props.VarDescription = catVarProps(a.props.VarDescription,b.props.VarDescription,LHSVars,RHSVars);
        a.props.Units = catVarProps(a.props.Units,b.props.Units,LHSVars,RHSVars);
    end

    if (maxObsIndex > a.nobs)
        % Don't warn if a had no variables originally
        if a.nvars > b.nvars
            warning(message('stats:dataset:subsasgn:DefaultValuesAdded'));
        end
        % Note that b's observation names do not propogate to a with ()
        % assignment.  They do for vertcat, though.
        a = fillInDataset(a,maxObsIndex,newObsNames);
    end

else
    % Raw values are not accepted as the RHS with '()' subscripting:  With a
    % single variable, you can use dot subscripting.  With multiple variables,
    % you can either wrap them up in a dataset, accepted above, or use braces
    % if the variables are homogeneous.
    error(message('stats:dataset:subsasgn:InvalidRHS'));
end
