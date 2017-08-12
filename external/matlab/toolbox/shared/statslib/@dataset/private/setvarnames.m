function a = setvarnames(a,newnames,vars,allowMods,allowDups,allowEmpty)
%SETVARNAMES Set dataset array variable names.

%   Copyright 2006-2012 The MathWorks, Inc.

if nargin < 4
    allowMods = true;
end
if nargin < 5
    allowDups = false;
end
if nargin < 6
    allowEmpty = false;
end

originalNames = newnames;

if isString(newnames)
    newnames = {newnames};
end

fullAssignment = (nargin == 2) || isempty(vars);
if fullAssignment % replacing all names
    varIndices = 1:a.nvars;
else % replacing some names
    varIndices = getvarindices(a,vars);
end

if ~isStrings(newnames,allowEmpty) % require a cell array, allow empty strings per caller
    error(message('stats:dataset:setvarnames:InvalidVarnames'));
elseif numel(newnames) ~= length(varIndices)
    if fullAssignment
        error(message('stats:dataset:setvarnames:IncorrectNumberOfVarnames'));
    else
        error(message('stats:dataset:setvarnames:IncorrectNumberOfVarnamesPartial'));
    end
end
newnames = strtrim(newnames(:))'; % a row vector, strtrim conveniently converts {} to a 1x0
[newnames,wasEmpty] = fixEmptyNames(newnames,varIndices);
[newnames, wasMadeValid] = matlab.lang.makeValidName(newnames);
if any(wasMadeValid) % will warn if mods are made
    if allowMods
        warning(message('stats:dataset:ModifiedVarnames'));
    else
        firstModified = newnames{find(wasMadeValid,1)};
        error(message('stats:dataset:InvalidVariableName',firstModified));
    end
end
if fullAssignment
    a_varnames = newnames;
else
    % Broadcast partial assignment out to full size
    a_varnames = a.varnames; a_varnames(varIndices) = newnames;
end
mods = false(size(a_varnames)); mods(varIndices(wasEmpty|wasMadeValid)) = true;

if allowDups
    % Uniqueify the new names (in their possibly modified form) with respect to
    % each other and to existing names
    a_varnames = matlab.lang.makeUniqueStrings(a_varnames,varIndices,namelengthmax);
else
    % Don't allow unmodified new names to duplicate existing names, but make
    % sure names that were filled in or made valid do not duplicate any other
    % names, either existing or new
    a_varnames = matlab.lang.makeUniqueStrings(a_varnames,mods,namelengthmax);
    checkduplicatenames(a_varnames,'varnames');
end
checkreservednames(a_varnames);
a.varnames = a_varnames;

if allowMods && any(wasMadeValid)
    vd = getproperty(a,'VarDescription',true);
    str = { getString(message('stats:dataset:uistrings:ModifiedVarnameDescr')) };
    vd(varIndices(wasMadeValid)) = strcat(str, {' '''}, originalNames(wasMadeValid), {''''});
    a = setvardescription(a,vd);
end


%-----------------------------------------------------------------------
function [names,empties] = fixEmptyNames(names,indices)
empties = cellfun('isempty',names);
if any(empties)
    names(empties) = dfltvarnames(indices(empties));
end


%-----------------------------------------------------------------------
function tf = isString(s)
% Require a (possibly empty) row of chars or ''.
tf = ischar(s) && (isrow(s) || isequal(s,''));


%-----------------------------------------------------------------------
function tf = isStrings(s,allowEmpty)
% ISSTRINGS Require a cell array of char row vectors, or ''
if allowEmpty
    stringTest = @(s) ischar(s) && ( isrow(s) || isequal(s,'') );
else
    stringTest = @(s) ischar(s) && isrow(s) && any(s ~= ' ');
end
if iscell(s)
    tf = all(cellfun(stringTest,s,'UniformOutput',true));
else
    tf = false;
end

