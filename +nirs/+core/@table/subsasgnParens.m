function t = subsasgnParens(t,s,b,creating)
%SUBSASGNPARENS Subscripted assignment to a table.

%   Copyright 2012-2014 The MathWorks, Inc.

import matlab.internal.tableUtils.matricize

% '()' is assignment to a subset of a table.  Only dot subscripting
% may follow.

if nargin < 4, creating = false; end
if ~isstruct(s), s = struct('type','()','subs',{s}); end

if numel(s(1).subs) ~= t.ndims
    error(message('MATLAB:table:NDSubscript'));
end

if ~isscalar(s)
    switch s(2).type
    case '()'
        error(message('MATLAB:table:InvalidSubscriptExpr'));
    case '{}'
        error(message('MATLAB:table:InvalidSubscriptExpr'));
    case '.'
        if creating
            error(message('MATLAB:table:InvalidSubscriptExpr'));
        end
        
        % Syntax:  t(rowIndices,varIndices).name = b
        % Syntax:  t(rowIndices,varIndices).name(...) = b
        % Syntax:  t(rowIndices,varIndices).name{...} = b
        % Syntax:  t(rowIndices,varIndices).name.field = b
        %
        % Assignment into a variable of a subarray.
        %
        % This may also be followed by deeper levels of subscripting.
        %
        % t(rowIndices,varIndices) must refer to rows and vars that exist, and
        % the .name assignment can't add rows or refer to a new variable.  This
        % is to prevent cases where the indexing beyond t(rowIndices,varIndices)
        % refers to things that are new relative to that subarray, but which
        % already exist in t itself.  So, cannot grow the table by an assignment
        % like this.
        %
        % This can be deletion, but it must be "inside" a variable, and not
        % change the size of t(rowIndices,varIndices).
        
        % Get the subarray, do the dot-variable assignment on that.
        try
            c = subsrefParens(t,s(1));
        catch ME
            outOfRangeIDs = {'MATLAB:table:RowIndexOutOfRange' 'MATLAB:table:UnrecognizedRowName' ...
                             'MATLAB:table:VarIndexOutOfRange' 'MATLAB:table:UnrecognizedVarName'};
            if any(strcmp(ME.identifier,outOfRangeIDs))
                error(message('MATLAB:table:InvalidExpansion'));
            else
                rethrow(ME);
            end
        end
        
        % Assigning to .Properties of a subarray is not allowed.
        if strcmp(s(2).subs,'Properties')
            error(message('MATLAB:table:PropertiesAssignmentToSubarray'));
        end
        
        nestedDeleting = builtin('_isEmptySqrBrktLiteral',b);
        b = subsasgnDot(c,s(2:end),b);
        
        % Changing the size of the subarray -- growing it by assignment or
        % deleting part of it -- is not allowed.
        if ~isequal(size(b),size(c))
            if nestedDeleting
                error(message('MATLAB:table:EmptyAssignmentToSubarrayVar'));
            else
                error(message('MATLAB:table:InvalidExpansion'));
            end
        end
        
        % Now let the simple () subscripting code handle assignment of the updated
        % subarray back into the original array.
        s = s(1);
    end
end

% If the RHS is (still) [], we are deleting a variable from the table.
deleting = builtin('_isEmptySqrBrktLiteral',b);

% If a new table is being created, or if the LHS is 0x0, then interpret
% ':' as the size of the corresponding dim from the RHS, not as nothing.
colonFromRHS = ~deleting && (creating || all(size(t)==0));

[b_nrows,b_nvars] = size(b);

% Translate row names into indices (leave ':' alone)
if colonFromRHS && iscolon(s(1).subs{1})
    rowIndices = 1:b_nrows;
    numRowIndices = b_nrows;
    maxRowIndex = b_nrows;
    newRowNames = {};
    isColonRows = true;
else
    [rowIndices,numRowIndices,maxRowIndex,newRowNames,isColonRows] = ...
                               getRowIndices(t, s(1).subs{1}, ~deleting);
end

% Translate variable (column) names into indices (translate ':' to 1:nvars)
if colonFromRHS && iscolon(s(1).subs{2})
    varIndices = 1:b_nvars;
    if iscell(b)
        newVarNames = matlab.internal.table.dfltVarNames(1:b_nvars);
    else
        newVarNames = b.varnames;
    end
else
    [varIndices,newVarNames] = getVarIndices(t, s(1).subs{2}, ~deleting);
end

% Syntax:  t(rowIndices,:) = []
%          t(:,varIndices) = []
%          t(rowIndices,varIndices) = [] is illegal
%
% Deletion of complete rows or entire variables.
if deleting
    % Delete rows across all variables
    if iscolon(s(1).subs{2})
        if isnumeric(rowIndices)
            rowIndices = unique(rowIndices);
            numRowIndices = numel(rowIndices);
        end
        newNrows = t.nrows - numRowIndices;
        t_data = t.data;
        for j = 1:t.nvars
            var_j = t_data{j};
            if isa(var_j,'table')
                var_j = subsasgnParens(var_j,{rowIndices ':'},[]); % can't use table subscripting directly
            elseif ismatrix(var_j)
                var_j(rowIndices,:) = []; % without using reshape, may not be one
            else
                sizeOut = size(var_j); sizeOut(1) = newNrows;
                var_j(rowIndices,:) = [];
                var_j = reshape(var_j,sizeOut);
            end
            t_data{j} = var_j;
        end
        t.data = t_data;
        if ~isempty(t.rownames), t.rownames(rowIndices) = []; end
        t.nrows = newNrows;

    % Delete entire variables
    elseif isColonRows
        varIndices = unique(varIndices); % getvarindices converts all varindex types to numeric
        t.data(varIndices) = [];
        t.varnames(varIndices) = [];
        t.nvars = t.nvars - numel(varIndices);
        % Var-based properties need to be shrunk.
        if ~isempty(t.props.VariableDescriptions), t.props.VariableDescriptions(varIndices) = []; end
        if ~isempty(t.props.VariableUnits), t.props.VariableUnits(varIndices) = []; end

    else
        error(message('MATLAB:table:InvalidEmptyAssignment'));
    end

% Syntax:  t(rowIndices,varIndices) = b
%
% Assignment from a table.  This operation is supposed to replace or
% grow at the level of the _table_.  So no internal reshaping of
% variables is allowed -- we strictly enforce sizes. In other words, the
% existing table has a specific size/shape for each variable, and
% assignment at this level must respect that.
else
    if isscalar(b) % scalar expansion of a single table element or cell, which may itself be non-scalar
        b = repmat(b,numRowIndices,length(varIndices));
        [b_nrows,b_nvars] = size(b);
    else
        if b_nrows ~= numRowIndices
            error(message('MATLAB:table:RowDimensionMismatch'));
        elseif b_nvars ~= length(varIndices)
            error(message('MATLAB:table:VarDimensionMismatch'));
        end
    end
    
    if isa(b,'table')
        b_data = b.data;
    elseif iscell(b)
        if ~ismatrix(b)
            error(message('MATLAB:table:NDCell'));
        end
        if b_nrows == 1
            b_data = b; % fast for row assignment
        else
            b_data = matlab.internal.table.container2vars(b);
        end
    else
        % Raw values are not accepted as the RHS with '()' subscripting:  With a
        % single variable, you can use dot subscripting.  With multiple variables,
        % you can either wrap them up in a table, accepted above, or use braces
        % if the variables are homogeneous.
        error(message('MATLAB:table:InvalidRHS'));
    end
    
    existingVars = (varIndices <= t.nvars);
    existingVarLocs = find(existingVars);
    t_data = t.data;
    [t_nrows,t_nvars] = size(t);
    for j = existingVarLocs
        var_j = t_data{varIndices(j)};
        % The size of the RHS has to match what it's going into.
        sizeLHS = size(var_j); sizeLHS(1) = numRowIndices;
        if prod(sizeLHS) ~= prod(size(b_data{j})) %#ok<PSIZE> avoid numel, it may return 1
            error(message('MATLAB:table:AssignmentDimensionMismatch', t.varnames{varIndices(j)}));
        end
        var_b = matricize(b_data{j});
% In cases where the whole var is moved, i.e. rowIndices is ':', this is faster, but a valid
% RHS may not have same type or trailing size as the LHS var, and it's difficult to do the
% right error checking - so do it as a subscripted assignment.
%         if isColonRows && isequal(sizeLHS,size(b_data{j}))) && isa(b_data{j},class(var_j))
%             var_j = var_b;
%         else
            if isa(var_j,'table')
                var_j = subsasgnParens(var_j,{rowIndices ':'},var_b); % can't use table subscripting directly
            else
                var_j(rowIndices,:) = var_b;
            end
%         end
        % No need to check for size change, RHS and LHS are identical sizes.
        t_data{varIndices(j)} = var_j;
    end

    % Add new variables if necessary.  Note that b's varnames do not
    % propagate to a in () assignment, unless t is being created or grown
    % from 0x0.  They do for horzcat, though.
    newVarLocs = find(~existingVars);
    if ~isempty(newVarLocs)
        matlab.internal.tableUtils.makeValidName(newVarNames,'error'); % error if any invalid
        
        % Warn if we have to lengthen the new variables to match the height of
        % the table. Don't warn about default values "filled in in the middle"
        % for these new vars.
        if maxRowIndex < t_nrows
            warning(message('MATLAB:table:RowsAddedNewVars'));
        end
        
        t_data = [t_data cell(1,length(newVarNames))];
        for j = 1:length(newVarNames)
            var_b = b_data{newVarLocs(j)};
            if isColonRows
                var_j = var_b;
            else
                % Start the new variable out as 0-by-(trailing size of b),
                % then let the assignment add rows.
                var_j = repmat(var_b,[0 ones(1,ndims(var_b)-1)]);
                if isa(var_j,'table')
                    var_j = subsasgnParens(var_j,{rowIndices ':'},matricize(var_b));
                else
                    var_j(rowIndices,:) = matricize(var_b);
                end
            end
            % A new var may need to grow to fit the table
            if size(var_j,1) < t_nrows
                var_j = lengthenVar(var_j, t_nrows);
            end
            t_data{t_nvars+j} = var_j;
        end
        t.varnames = [t.varnames newVarNames];
        t.nvars = t_nvars + length(newVarNames);
        % Var-based properties need to be extended.
        LHSVars = 1:t_nvars;
        RHSVars = 1:b_nvars;
        if iscell(b)
            t.props.VariableDescriptions = catVarProps(t.props.VariableDescriptions,{},LHSVars,RHSVars);
            t.props.VariableUnits = catVarProps(t.props.VariableUnits,{},LHSVars,RHSVars);
        else
            t.props.VariableDescriptions = catVarProps(t.props.VariableDescriptions,b.props.VariableDescriptions,LHSVars,RHSVars);
            t.props.VariableUnits = catVarProps(t.props.VariableUnits,b.props.VariableUnits,LHSVars,RHSVars);
        end
    end
    t.data = t_data;

    if (maxRowIndex > t_nrows) % t_nrows is t's original number of rows
        % If the vars being assigned to are now taller than the table, add rows
        % to the rest of the table, including row names.  This might be because
        % the assignment lengthened existing vars, or because the assignment
        % created new vars taller than the table.  Warn only if we have to
        % lengthen existing vars that have not been assigned to -- if there's
        % currently only one var in the table (which might be existing or new),
        % don't warn about any default values "filled in in the middle".
        if length(varIndices) < t.nvars; % some existing vars were not assigned to
            warning(message('MATLAB:table:RowsAddedExistingVars'));
        end
        % Note that b's row names do not propagate to t with ()
        % assignment.  They do for vertcat, though.
        t = fillOutTable(t,maxRowIndex,newRowNames);
    end
end
