function [Fpass1, Fstop1, Fstop2, Fpass2, d1, d2, d3] = getdesignspecs(h, d);
%GETDESIGNSPECS

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Get frequency specs, they have been prenormalized
Fpass1 = get(d,'Fpass1');
Fstop1 = get(d,'Fstop1');
Fstop2 = get(d,'Fstop2');
Fpass2 = get(d,'Fpass2');

% Set the magUnits temporarily to 'linear' to get deviations
magUnits = get(d,'magUnits');
set(d,'magUnits','linear');
d1 = get(d,'Dpass1');
d2 = get(d,'Dstop');
d3 = get(d,'Dpass2');
set(d,'magUnits',magUnits);

% [EOF]
