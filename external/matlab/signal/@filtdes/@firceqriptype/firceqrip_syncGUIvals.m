function firceqrip_syncGUIvals(h, d, arrayh)
%FIRCEQRIP_SYNCGUIVALS Sync values from frames.
%
%   Inputs:
%       h - handle to this object   
%       arrayh - array of handles to frames

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Get handles to options frame

ft_syncGUIvals(h, d, arrayh);

fr    = getoptsframe(h);
hopts = find(arrayh,'-class', fr.constructor);

% Store specs in object
set(d,'minPhase',get(hopts,'isMinPhase'));
set(d,'stopbandSlope',evaluatevars(get(hopts,'StopbandSlope')));

% [EOF]
