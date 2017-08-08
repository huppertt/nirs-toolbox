function [Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2] = getdesignspecs(h,d);
%GETDESIGNSPECS

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Get frequency specs, they have been prenormalized
Fstop1 = get(d,'Fstop1');
Fpass1 = get(d,'Fpass1');
Fpass2 = get(d,'Fpass2');
Fstop2 = get(d,'Fstop2');

% Set the magUnits temporarily to 'dB' to get attenuations
magUnits = get(d,'magUnits');
set(d,'magUnits','dB');
Astop1 = get(d,'Astop1');
Apass = get(d,'Apass');
Astop2 = get(d,'Astop2');
set(d,'magUnits',magUnits);

% [EOF]
