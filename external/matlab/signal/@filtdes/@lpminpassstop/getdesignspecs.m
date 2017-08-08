function [Fpass, Fstop, delta1, delta2] = getdesignspecs(hObj, d)
%GETDESIGNSPECS

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Get frequency specs, they have been prenormalized
Fpass = get(d,'Fpass');
Fstop = get(d,'Fstop');

% Set the magUnits temporarily to 'linear' to get deviations
magUnits = get(d,'magUnits');
set(d,'magUnits','linear');
delta1 = get(d,'Dpass');
delta2 = get(d,'Dstop');
set(d,'magUnits',magUnits);

% [EOF]
