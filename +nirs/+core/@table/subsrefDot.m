function [varargout] = subsrefDot(t,s)
%SUBSREFDOT Subscripted reference for a table.

%   Copyright 2013-2014 The MathWorks, Inc.

import matlab.internal.tableUtils.isstring

% '.' is a reference to a table variable or property.  Any sort of
% subscripting may follow.  Row names for cascaded () or {} subscripting on
% a variable are inherited from the table.

if ~isstruct(s), s = struct('type','.','subs',s); end

% Translate variable (column) name into an index.
varName = s(1).subs;
if isnumeric(varName) && isscalar(varName)
    % Allow t.(i) where i is an integer
    varIndex = varName;
elseif isstring(varName)
    varIndex = find(strcmp(varName,t.varnames));
    if isempty(varIndex)
        % If there's no such var, it may be a reference to the 'properties'
        % (virtual) property.  Handle those, but disallow references to
        % any property directly.
        if strcmp(varName,'Properties')
            if isscalar(s)
                varargout{1} = getProperties(t);
            else
                % If there's cascaded subscripting into the property, let the
                % property's subsref handle the reference. This may result in
                % a comma-separated list, so ask for and assign to as many
                % outputs as we're given. That is the number of outputs on
                % the LHS of the original expression, or if there was no LHS,
                % it comes from numArgumentsFromSubscript.
                try
                    [varargout{1:nargout}] = getProperty(t,s(2:end));
                catch ME
                    if strcmp(ME.identifier,'MATLAB:table:UnknownProperty')
                        propName = s(2).subs;
                        match = find(strcmpi(propName,table.propertyNames),1);
                        if ~isempty(match) % a property name, but with wrong case
                            match = table.propertyNames{match};
                            error(message('MATLAB:table:UnknownPropertyCase',propName,match));
                        else
                            throw(ME);
                        end
                    else
                        throw(ME);
                    end
                end
            end
            return
        elseif strcmpi(varName,'Properties') % .Properties, but with wrong case
            error(message('MATLAB:table:UnrecognizedVarNamePropertiesCase',varName));
        else
            match = find(strcmpi(varName,table.propertyNames),1);
            if ~isempty(match)
                match = table.propertyNames{match};
                if strcmp(varName,match) % a valid property name
                    error(message('MATLAB:table:IllegalPropertyReference',varName,varName));
                else % a property name, but with wrong case
                    error(message('MATLAB:table:IllegalPropertyReferenceCase',varName,match,match));
                end
            else
                match = find(strcmpi(varName,t.varnames),1);
                if ~isempty(match) % an existing variable name
                    match = t.varnames{match};
                    error(message('MATLAB:table:UnrecognizedVarNameCase',varName,match));
                else
                    methodList = methods(t);
                    match = find(strcmpi(varName,methodList),1);
                    if ~isempty(match) % a method name
                        match = methodList{match};
                        error(message('MATLAB:table:IllegalDotMethod',varName,match,match));
                    else % no obvious match
                        error(message('MATLAB:table:UnrecognizedVarName',varName));
                    end
                end
            end
        end
    end
else
    error(message('MATLAB:table:IllegalVarSubscript'));
end

b = t.data{varIndex};

if isscalar(s)
    % If there's no additional subscripting, return the table variable.
    varargout{1} = b;
else
    if ~isequal(s(2).type,'.') % () or {} subscripting after dot
        % The variable inherits row names from the table.
        % Translate names to row numbers if necessary.
        rowIndices = s(2).subs{1};
        if iscolon(rowIndices) || islogical(rowIndices) || isnumeric(rowIndices)
            % leave these alone
        else
            numericRowIndices = getRowIndices(t, rowIndices); % (leaves ':' alone)
            if (size(b,2)>1) && isscalar(s(2).subs)
                error(message('MATLAB:table:InvalidLinearIndexing'));
            end
            % getRowIndices returns the indices as a col vector, but subscripting on
            % a table variable (as opposed to on a table) should follow the usual
            % reshaping rules. Nothing to do for one (char) name, but preserve a
            % cellstr subscript's original shape.
            if iscell(rowIndices), numericRowIndices = reshape(numericRowIndices,size(rowIndices)); end
            s(2).subs{1} = numericRowIndices;
        end
    else
        % A reference to a property or field, so no row labels
    end
    
    % Now let the variable's subsref handle the remaining subscripts in things
    % like t.name(...) or  t.name{...} or t.name.property. This may return a
    % comma-separated list, so ask for and assign to as many outputs as we're
    % given. That is the number of outputs on the LHS of the original expression,
    % or if there was no LHS, it comes from numArgumentsFromSubscript.
    if length(s) == 2
        try %#ok<ALIGN>
            [varargout{1:nargout}] = subsref(b,s(2));
        catch ME, throw(ME); end
    else % length(s) > 2
        % Trick the third and higher levels of subscripting in things like
        % t.Var{i}(...) etc. into dispatching to the right place when
        % t.Var{i}, or something further down the chain, is itself a table.
        try %#ok<ALIGN>
            [varargout{1:nargout}] = matlab.internal.table.subsrefRecurser(b,s(2:end));
        catch ME, rethrow(ME); end % point to the line in subsrefRecurser
    end
end
