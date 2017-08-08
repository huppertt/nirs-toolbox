function [Fstop1, Fstop2, Astop] = getdesignspecs(h, d)
%GETDESIGNSPECS Return the specs for the design

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Get frequency specs, they have been prenormalized
Fstop1 = get(d,'Fstop1');
Fstop2 = get(d,'Fstop2');

% Set the magUnits temporarily to 'dB' to get attenuation
magUnits = get(d,'magUnits');
set(d,'magUnits','dB');
Astop = get(d,'Astop');
set(d,'magUnits',magUnits);

% [EOF]
