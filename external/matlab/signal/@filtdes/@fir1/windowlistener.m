function windowlistener(this,eventdata)
%WINDOWLISTENER Callback for listener to the window property.
    
%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

rmprops(this, this.WindowObject);

% Get window
win = get(this,'Window');

winobj = getwinobject(win);

% If window is kaiser and minimum order is supported for this type,
% enable ordermode property otherwise disable it                   
if isa(winobj, 'sigwin.kaiser') & ~isempty(findConstr(this,get(this,'responseType'),'minimum')),
    enab = 'on';
else
    enab = 'off';
end
enabdynprop(this,'orderMode',enab);

set(this, 'WindowObject', winobj);

prop = getwindowprop(this);
if ~isempty(prop),
    addprops(this, this.WindowObject, prop);
end

% [EOF]
