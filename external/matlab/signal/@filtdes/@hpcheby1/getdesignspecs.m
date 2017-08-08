function [Fpass, Apass] = getdesignspecs(h, d)
%GETDESIGNSPECS Return the specs for design

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Get frequency specs, they have been prenormalized
Fpass = get(d,'Fpass');

% Set the magUnits temporarily to 'dB' to get attenuation
magUnits = get(d,'magUnits');
set(d,'magUnits','dB');
Apass = get(d,'Apass');
set(d,'magUnits',magUnits);

% [EOF]
