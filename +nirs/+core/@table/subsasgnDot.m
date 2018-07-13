function t = subsasgnDot(t,s,b)
%SUBSASGNDOT Subscripted assignment to a table.

%   Copyright 2013-2014 The MathWorks, Inc.

% '.' is assignment to or into a variable.  Any sort of subscripting
% may follow that, and row names are inherited from the table.

import matlab.internal.tableUtils.emptyLike
import matlab.internal.tableUtils.isstring

if ~isstruct(s), s = struct('type','.','subs',s); end

% For built-in tyes, s(2).type=='()' guarantees that length(s)==2, but for a var that
% is itself a table, parens need not be the end. deletion other than t.Var or t.Var(...)
% is handled by the assignment code path.
deleting = builtin('_isEmptySqrBrktLiteral',b) ...
    && (isscalar(s) || ((length(s) == 2) && isequal(s(2).type,'()')));

% Translate variable (column) name into an index.
varName = s(1).subs;
if isnumeric(varName) && isscalar(varName)
    % Allow t.(i) where i is an integer
    [varIndex,newNames] = getVarIndices(t,varName,~deleting);
    isNewVar = (varIndex > t.nvars);
    if isNewVar
        varName = newNames{1};
    else
        varName = t.varnames{varIndex};
    end
elseif isstring(varName)
    varIndex = find(strcmp(varName,t.varnames));
    isNewVar = isempty(varIndex);
    if isNewVar
        if deleting
            error(message('MATLAB:table:UnrecognizedVarName', varName));
        end
        % Handle assignment to a property under the 'Properties' (virtual)
        % property, or to the entire 'Properties' property.
        if matlab.internal.table.checkReservedNames(varName)
            if strcmp(varName,'Properties')
                try
                    if isscalar(s)
                        t = setProperties(t,b);
                    else
                        t = setProperty(t,s(2:end),b);
                    end
                catch ME
                    if ~isscalar(s) && strcmp(ME.identifier,'MATLAB:table:UnknownProperty')
                        propName = s(2).subs;
                        match = find(strcmpi(propName,table.propertyNames),1);
                        if ~isempty(match) % a property name. but with wrong case
                            match = table.propertyNames{match};
                            error(message('MATLAB:table:UnknownPropertyCase',propName,match));
                        else
                            throw(ME);
                        end
                    else
                        throw(ME);
                    end
                end
                return
            else % t.VariableNames or t.RowNames
                error(message('MATLAB:table:InvalidPropertyAssignment',varName,varName));
            end
        end
        matlab.internal.tableUtils.makeValidName(varName,'error'); % error if invalid
        
        % If this is a new variable, it will go at the end.
        varIndex = t.nvars + 1;
    end
else
    error(message('MATLAB:table:IllegalVarSubscript'));
end

% Handle empty assignment intended as deletion of an entire variable or of
% columns/pages/etc. of a variable.  Deletion of rows in a (single)
% variable is caught here and not allowed.  Other empty assignment
% syntaxes may be assignment to cells or may be deletion of things deeper
% in a non-atomic variable, neither is handled here.
if deleting
    % Syntax:  t.var = []
    %
    % Delete an entire variable.
    if isscalar(s)
        t.data(varIndex) = [];
        t.varnames(varIndex) = [];
        t.nvars = t.nvars - 1;
        % Var-based or properties need to be shrunk.
        if ~isempty(t.props.VariableDescriptions), t.props.VariableDescriptions(varIndex) = []; end
        if ~isempty(t.props.VariableUnits), t.props.VariableUnits(varIndex) = []; end
        
    % Syntax:  t.var(:,...) = []
    %          t.var(rowIndices,...) = [] is illegal
    %
    % Delete columns/pages/etc. of a variable, with ':' as the first index
    % in subscript.  This may change the dimensionality of the variable,
    % but won't change the number of rows because we require ':' as the
    % first index.
    else % length(s) == 2
        if ~iscolon(s(2).subs{1})
            error(message('MATLAB:table:EmptyAssignmentToVariableRows'));
        elseif isscalar(s(2).subs)
            error(message('MATLAB:table:EmptyAssignmentToAllVariableElements'));
        end
        
        var_j = t.data{varIndex};
        try %#ok<ALIGN>
            if isa(var_j,'table')
                var_j = subsasgnParens(var_j,s(2),[]); % can't use table subscripting directly
            else
                var_j(s(2).subs{:}) = [];
            end
        catch ME, throw(ME); end
        t.data{varIndex} = var_j;
    end
    
else
    % Syntax:  t.var = b
    %
    % Replace an entire variable.  It must have the right number of rows, unless
    % the RHS is 0x0.
    if isscalar(s)
        if size(b,1) ~= t.nrows && ~all(size(t) == 0)
            % If the assignment has the wrong number of rows, check for some
            % common mistakes to suggest what may have been intended
            if strcmpi(varName,'Properties') && isstruct(b) && isscalar(b)
                % Things like t.properties = scalarStruct
                str = getString(message('MATLAB:table:IntendedPropertiesAssignment'));
                error(message('MATLAB:table:RowDimensionMismatchSuggest',str));
            else
                match = find(strcmpi(varName,table.propertyNames),1);
                if ~isempty(match)
                    % Things like t.PropertyName = ...
                    match = table.propertyNames{match};
                    str = getString(message('MATLAB:table:IntendedPropertyAssignment',match,match));
                    error(message('MATLAB:table:RowDimensionMismatchSuggest',str));
                end
            end
            % Anything else, no suggestion. No point in checking for a case
            % insensitive match to an existing var, even with the correct case,
            % this would still be an illegal assignment
            error(message('MATLAB:table:RowDimensionMismatch'));
        end
        var_j = b;
        newRowNames = {};
        
    % Syntax:  t.var(rowIndices,...) = b
    %          t.var{rowIndices,...} = b
    %          t.var{rowIndices,...} = [] (this is assignment, not deletion)
    %          t.var.field = b
    %
    % Assign to elements in a variable.  Assignment can also be used to
    % expand the variable's number of rows, or along another dimension.
    %
    % Cell indexing, e.g. t.var{rowIndices,...}, or a reference to a
    % field, e.g. t.var.field, may also be followed by deeper levels of
    % subscripting. Cannot create a new var implicitly by deeper indexing.
    else % ~isscalar(s)
        if isNewVar && (length(s) > 2) && ~isequal(s(2).type,'.')
            % If the assignment is not to an existing var, check for some common
            % mistakes to suggest what may have been intended
            match = find(strcmpi(varName,t.varnames),1);
            if ~isempty(match)
                % An existing variable name, but with wrong case
                match = t.varnames{match};
                str = getString(message('MATLAB:table:IntendedVarAssignment',match));
                error(message('MATLAB:table:InvalidExpansionDotDepthSuggest',str));
            end
            % Anything else, no suggestion
            error(message('MATLAB:table:InvalidExpansionDotDepth'));
        end
        if isequal(s(2).type,'.') % dot indexing into variable
            % If the assignment is not to an existing var, check for some common
            % mistakes to suggest what may have been intended
            if isNewVar
                if strcmpi(varName,'Properties') && isstring(s(2).subs)
                    % Things like t.properties.name
                    str = getString(message('MATLAB:table:IntendedPropertiesAssignment'));
                    error(message('MATLAB:table:InvalidExpansionDotSuggest',str));
                else
                    match = find(strcmpi(varName,t.varnames),1);
                    if ~isempty(match)
                        % An existing variable name, but with wrong case
                        match = t.varnames{match};
                        str = getString(message('MATLAB:table:IntendedVarAssignment',match));
                        error(message('MATLAB:table:InvalidExpansionDotSuggest',str));
                    else
                        % Anything else, no suggestion
                        error(message('MATLAB:table:InvalidExpansionDot'));
                    end
                end
            end
            var_j = t.data{varIndex};
        else % () or {} subscripting into variable
            % Initialize a new var, or extract an existing var.
            if isNewVar
                % Start the new var out as an empty of b's class with
                % the same number of rows as the table.
                var_j = emptyLike([t.nrows,0],'Like',b);
            else
                var_j = t.data{varIndex};
            end
            
            % The variable inherits row names from the table.  Translate names to
            % row numbers if necessary.
            rowIndices = s(2).subs{1};
            if iscolon(rowIndices) || islogical(rowIndices) || isnumeric(rowIndices)
                % The variable can whatever it wants with non-rowName indices.
            else
                % Row names may be used as a "linear index" only if the variable
                % has a single column.
                if ~iscolumn(var_j) && isscalar(s(2).subs)
                    error(message('MATLAB:table:InvalidLinearIndexing'));
                end
            end
            % getRowIndices returns the indices as a col vector, which prevents reshaping.
            % This is fine because the var is constrained inside the table. 
            [rowIndices,~,~,newRowNames] = getRowIndices(t, rowIndices, true);
            s(2).subs{1} = rowIndices;
        end
        
        % Now let the variable's subsasgn handle the subscripting in
        % things like t.name(...) or  t.name{...} or t.name.attribute
        
        % Calling subsasgn directly allows changing the shape of var_j by
        % assignment, which the interpreter would not allow.
        if length(s) == 2
            if isobject(var_j)
                % If b is superior to var_j, this will dispatch to b
                try %#ok<ALIGN>
                    var_j = subsasgn(var_j,s(2),b);
                catch ME, throw(ME); end
            else
                % Call builtin, to get correct dispatching even if b is an object.
                try %#ok<ALIGN>
                    var_j = builtin('subsasgn',var_j,s(2),b);
                catch ME, throw(ME); end
            end
        else % length(s) > 2
            % Trick the third and higher levels of subscripting in things like
            % t.Var{i}(...) etc. into dispatching to the right place even when
            % t.Var{i}, or something further down the chain, is itself a table.
            try %#ok<ALIGN>
                var_j = matlab.internal.table.subsasgnRecurser(var_j,s(2:end),b);
            catch ME, rethrow(ME); end % point to the line in subsasgnRecurser
        end
    end
    
    % If this is a new variable, make it official.
    if isNewVar
        t.varnames = [t.varnames varName];
        t.nvars = varIndex;
        % Var-based properties need to be extended.
        if ~isempty(t.props.VariableDescriptions), t.props.VariableDescriptions = [t.props.VariableDescriptions {''}]; end
        if ~isempty(t.props.VariableUnits), t.props.VariableUnits = [t.props.VariableUnits {''}]; end
    end
    
    % If an entire var was replaced or created, the new value was required to
    % have the same number of rows as the table.  However, when assigning into a
    % new var, the assignment might create something shorter than the table, so
    % check for that and lengthen the new var to match the table.  In most
    % cases, assigning into an existing var leaves it as the correct height,
    % although in at least one oddball case (assigning a field to a non-struct),
    % assigning into an existing var can make the result shorter then it was
    % originally, so check for those too.
    varLen = size(var_j,1);
    if varLen < t.nrows
        warning(message('MATLAB:table:RowsAddedNewVars'));
        var_j = lengthenVar(var_j,t.nrows);
    end
    t.data{varIndex} = var_j;
    
    % If the var being assigned to is now taller than the table, add rows to
    % the rest of the table, including row names.  This might be because the
    % assignment lengthened an existing var, or because an "into" assignment
    % created a new var taller than the table.  Warn only if we have to lengthen
    % existing vars that have not been assigned to -- if there's currently only
    % one var in the table (which might be existing or new), don't warn about
    % any default values "filled in in the middle".
    if varLen > t.nrows % t.nrows still the original value, has not grown yet
        if 1 < t.nvars; % some existing vars were not assigned to
            warning(message('MATLAB:table:RowsAddedExistingVars'));
        end
        t = fillOutTable(t,varLen,newRowNames); % updates nrows
    end
end
