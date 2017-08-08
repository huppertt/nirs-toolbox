function [Fpass1, Fpass2, Apass] = getdesignspecs(h, d)
%GETDESIGNSPECS Returns the specs for the design

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Get frequency specs, they have been prenormalized
Fpass1 = get(d,'Fpass1');
Fpass2 = get(d,'Fpass2');

% Set the magUnits temporarily to 'dB' to get attenuation
magUnits = get(d,'magUnits');
set(d,'magUnits','dB');
Apass = get(d,'Apass');
set(d,'magUnits',magUnits);

% [EOF]
