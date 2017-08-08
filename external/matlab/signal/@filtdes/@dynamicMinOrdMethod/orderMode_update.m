function orderMode_update(h,ft)
%ORDERMODE_UPDATE Update the order mode property.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

constr = findConstr(h,ft,'minimum');

if isempty(constr),
    % Disable orderMode
    enabdynprop(h,'orderMode','off');
else
    % Enable orderMode
    enabdynprop(h,'orderMode','on');
end
