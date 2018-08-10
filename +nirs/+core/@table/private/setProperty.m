function t = setProperty(t,name,p)
%SETPROPERTY Set a table property.

%   Copyright 2012 The MathWorks, Inc.

import matlab.internal.tableUtils.isStrings

% We may be given a name (when called from set), or a subscript expression
% that starts with a '.name' subscript (when called from subsasgn).  Get the
% name and validate it in any case.
if isstruct(name)
    s = name;
    if s(1).type ~= '.'
        error(message('MATLAB:table:InvalidSubscript'));
    end
    name = s(1).subs;
    haveSubscript = true;
else
    haveSubscript = false;
end
% Allow partial match for property names if this is via the set method;
% require exact match if it is direct assignment via subsasgn
name = matchPropertyName(name,haveSubscript);

if haveSubscript && ~isscalar(s)
    % If This is 1-D named parens/braces subscripting, convert the names to 
    % correct indices for properties that support named indexing. 
    % e.g. t.Properties.RowNames('SomeRowName')
    if ~strcmp(s(2).type,'.') && isscalar(s(2).subs)
        sub = s(2).subs{1};
        if isStrings(sub) % a name, names, or colon
            switch name
            case {'VariableNames' 'VariableDescriptions' 'VariableUnits'}
                s(2).subs{1} = getVarIndices(t,sub);
            case 'RowNames'
                s(2).subs{1} = getRowIndices(t,sub);
            case 'DimensionNames'
                s(2).subs{1} = getDimIndices(t,sub);
            end
        end
    end
    % If there's cascaded subscripting into the property, get the existing
    % property value and let the property's subsasgn handle the assignment.
    % This may change its shape or size or otherwise make it invalid; that
    % gets checked by the individual setproperty methods called below.  The
    % property may currently be empty, ask for a non-empty default version to
    % allow assignment into only some elements.
    oldp = getProperty(t,name,true);
    p = subsasgn(oldp,s(2:end),p);
end

% Assign the new property value into the dataset.
switch name
case 'RowNames'
    t = setRowNames(t,p);
case 'VariableNames'
    t = setVarNames(t,p); % error if invalid, duplicate, or empty
case 'DimensionNames'
    t = setDimNames(t,p);
case 'VariableDescriptions'
    t = setVarDescription(t,p);
case 'VariableUnits'
    t = setUnits(t,p);
case 'Description'
    t = setDescription(t,p);
case 'UserData'
    t = setUserData(t,p);
end
