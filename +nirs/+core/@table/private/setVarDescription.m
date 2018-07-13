function t = setVarDescription(t,newvardescr)
%SETVARDESCRIPTION Set table VariableDescriptions property.

%   Copyright 2012 The MathWorks, Inc.

import matlab.internal.tableUtils.isStrings

if ~isStrings(newvardescr,true) % require a cell array, allow empty strings
    error(message('MATLAB:table:InvalidVarDescr'));
elseif ~isempty(newvardescr) && numel(newvardescr) ~= t.nvars
    error(message('MATLAB:table:IncorrectNumberOfVarDescrs'));
end

if isequal(size(newvardescr),[1 0]) && t.nvars == 0
    % leave a 1x0 cell alone for a table with no vars
elseif isempty(newvardescr)
    newvardescr = {}; % for cosmetics
else
    newvardescr = strtrim(newvardescr(:))'; % a row vector
end
t.props.VariableDescriptions = newvardescr;
