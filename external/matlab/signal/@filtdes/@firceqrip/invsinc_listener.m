function invsinc_listener(d,eventdata)
%INVSINC_LISTENER Callback for listener to the invSinc property.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Determine if invSinc is on or off
invsincstate = get(d,'invSinc');

% Enable/disable the relevant dynamic properties
enabdynprop(d,'invSincFreqFactor',invsincstate);
enabdynprop(d,'invSincPower',invsincstate);