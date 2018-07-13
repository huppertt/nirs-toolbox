function t = setUnits(t,newunits)
%SETUNITS Set table VariableUnits property.

%   Copyright 2012 The MathWorks, Inc.

import matlab.internal.tableUtils.isStrings

if ~isStrings(newunits,true) % require a cell array, allow empty strings
    error(message('MATLAB:table:InvalidUnits'));
elseif ~isempty(newunits) && numel(newunits) ~= t.nvars
    error(message('MATLAB:table:IncorrectNumberOfUnits'));
end

if isequal(size(newunits),[1 0]) && t.nvars == 0
    % leave a 1x0 cell alone for a table with no vars
elseif isempty(newunits)
    newunits = {}; % for cosmetics
else
    newunits = strtrim(newunits(:))'; % a row vector
end
t.props.VariableUnits = newunits;
