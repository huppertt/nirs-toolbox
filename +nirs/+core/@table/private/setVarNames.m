function t = setVarNames(t,newnames,vars,allowMods,allowDups,allowEmpty)
%SETVARNAMES Set table variable names.

%   Copyright 2012 The MathWorks, Inc.

import matlab.internal.tableUtils.isStrings
import matlab.internal.tableUtils.isstring

if nargin < 4
    allowMods = false;
end
if nargin < 5
    allowDups = false;
end
if nargin < 6
    allowEmpty = false;
end

originalNames = newnames;

fullAssignment = (nargin == 2) || builtin('_isEmptySqrBrktLiteral',vars);
if fullAssignment % replacing all names
    varIndices = 1:t.nvars;
else % replacing some names
    if isstring(newnames)
        newnames = {newnames};
    end
    varIndices = getVarIndices(t,vars);
end

if ~isStrings(newnames,true,allowEmpty) % require a cell array, allow empty strings per caller
    error(message('MATLAB:table:InvalidVarNames'));
elseif numel(newnames) ~= length(varIndices)
    if fullAssignment
        error(message('MATLAB:table:IncorrectNumberOfVarNames'));
    else
        error(message('MATLAB:table:IncorrectNumberOfVarNamesPartial'));
    end
end
newnames = strtrim(newnames(:))'; % a row vector, strtrim conveniently converts {} to a 1x0
[newnames,wasEmpty] = matlab.internal.table.fixEmptyNames(newnames,varIndices);

if allowMods
    exceptionMode = 'warn';
else
    exceptionMode = 'error';
end
[newnames,wasMadeValid] = matlab.internal.tableUtils.makeValidName(newnames,exceptionMode); % will warn if mods are made

if fullAssignment
    t_varnames = newnames;
else
    % Broadcast partial assignment out to full size
    t_varnames = t.varnames; t_varnames(varIndices) = newnames;
end
mods = false(size(t_varnames)); mods(varIndices(wasEmpty|wasMadeValid)) = true;

if allowDups
    % Uniqueify the new names (in their possibly modified form) with respect to
    % each other and to existing names
    t_varnames = matlab.lang.makeUniqueStrings(t_varnames,varIndices,namelengthmax);
else
    % Don't allow unmodified new names to duplicate existing names, but make
    % sure names that were filled in or made valid do not duplicate any other
    % names, either existing or new
    t_varnames = matlab.lang.makeUniqueStrings(t_varnames,mods,namelengthmax);
    checkDuplicateNames(t_varnames,'varnames');
end
matlab.internal.table.checkReservedNames(t_varnames);
t.varnames = t_varnames;

if allowMods && any(wasMadeValid)
    vd = getProperty(t,'VariableDescriptions',true);
    str = { getString(message('MATLAB:table:uistrings:ModifiedVarNameDescr')) };
    vd(varIndices(wasMadeValid)) = strcat(str, {' '''}, originalNames(wasMadeValid), {''''});
    t = setVarDescription(t,vd);
end