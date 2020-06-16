function t = setDimNames(t,newnames,allowDups,allowEmpty)
%SETDIMNAMES Set table DimensionNames property.

%   Copyright 2012 The MathWorks, Inc.

import matlab.internal.tableUtils.isStrings

if nargin < 3
    allowDups = false;
end
if nargin < 4
    allowEmpty = false;
end

if ~isStrings(newnames,true,allowEmpty) % require a cell array, allow empty strings per caller
    error(message('MATLAB:table:InvalidDimNames'));
elseif numel(newnames) ~= t.ndims
    error(message('MATLAB:table:IncorrectNumberOfDimNames'));
end
[t_dimnames,mods] = fixEmptyNames(newnames);

if allowDups
    t_dimnames = matlab.lang.makeUniqueStrings(t_dimnames,{},namelengthmax);
else
    % Don't allow duplicate names, but make sure empty names that were filled in
    % do not duplicate any other names
    t_dimnames = matlab.lang.makeUniqueStrings(t_dimnames,mods,namelengthmax);
    checkDuplicateNames(t_dimnames,'dimnames');
end
t.props.DimensionNames = t_dimnames(:)'; % a row vector


%-----------------------------------------------------------------------
function [names,empties] = fixEmptyNames(names)
empties = cellfun('isempty',names);
if any(empties)
    dfltNames = matlab.internal.table.dfltDimNames();
    names(empties) = dfltNames(empties);
end
