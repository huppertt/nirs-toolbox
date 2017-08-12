function a = setunits(a,newunits)
%SETUNITS Set dataset array Units property.

%   Copyright 2006-2011 The MathWorks, Inc.


if nargin < 2
    error(message('stats:dataset:setunits:TooFewInputs'));
end

if nargin == 2
    if isempty(newunits)
        a.props.Units = {}; % do this for cosmetics
        return
    end
    if numel(newunits) ~= a.nvars
        error(message('stats:dataset:setunits:WrongLength'));
    elseif ~iscell(newunits)
        error(message('stats:dataset:setunits:InvalidUnits'));
    elseif ~all(cellfun(@isstring,newunits))
        error(message('stats:dataset:setunits:InvalidUnits'));
    end
    a.props.Units = newunits(:)';
end

function tf = isstring(s) % require a row of chars, or possibly ''
tf = ischar(s) && ((isvector(s) && (size(s,1) == 1)) || all(size(s)==0));
