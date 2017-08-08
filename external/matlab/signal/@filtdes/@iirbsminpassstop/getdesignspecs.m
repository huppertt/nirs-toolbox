function [Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, Apass2] = getdesignspecs(h,d)
%GETDESIGNSPECS

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Get frequency specs, they have been prenormalized
Fpass1 = get(d,'Fpass1');
Fstop1 = get(d,'Fstop1');
Fstop2 = get(d,'Fstop2');
Fpass2 = get(d,'Fpass2');

% Set the magUnits temporarily to 'dB' to get attenuations
magUnits = get(d,'magUnits');
set(d,'magUnits','dB');
Apass1 = get(d,'Apass1');
Astop = get(d,'Astop');
Apass2 = get(d,'Apass2');
set(d,'magUnits',magUnits);
