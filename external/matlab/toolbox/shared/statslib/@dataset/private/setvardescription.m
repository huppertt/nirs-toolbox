function a = setvardescription(a,newvardescr)
%SETVARDESCRIPTION Set dataset array VarDescription property.

%   Copyright 2006-2011 The MathWorks, Inc.


if nargin < 2
    error(message('stats:dataset:setvardescription:TooFewInputs'));
end

if nargin == 2
    if isempty(newvardescr)
        a.props.VarDescription = {}; % do this for cosmetics
        return
    end
    if numel(newvardescr) ~= a.nvars
        error(message('stats:dataset:setvardescription:WrongLength'));
    elseif ~iscell(newvardescr)
        error(message('stats:dataset:setvardescription:InvalidDescr'));
    elseif ~all(cellfun(@isstring,newvardescr))
        error(message('stats:dataset:setvardescription:InvalidDescr'));
    end
    a.props.VarDescription = newvardescr(:)';
end

function tf = isstring(s) % require a row of chars, or possibly ''
tf = ischar(s) && ((isvector(s) && (size(s,1) == 1)) || all(size(s)==0));
