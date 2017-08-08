function d = fir1
%FIRLS  Constructor for this design method object.
%
%   Outputs:
%       d - Handle to the design method object

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

d = filtdes.fir1;

% Call super's constructor
dynMinOrder_construct(d);

% Setup a listener for the window property
l = handle.listener(d,findprop(d,'Window'),'PropertyPostSet',@windowlistener);
% Store listener
set(l, 'callbacktarget', d); % Allow for methods as callbacks
set(d,'windowlistener',l);

% Fire the listeners manually the first time
windowlistener(d);
update_winparam(d);

% Set the tag
set(d,'Tag','FIR Window');

