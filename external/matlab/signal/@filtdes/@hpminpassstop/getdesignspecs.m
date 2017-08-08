function [Fstop, Fpass, delta1, delta2] = getdesignspecs(hObj, d)
%GETDESIGNSPECS

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

Fstop = get(d,'Fstop');
Fpass = get(d,'Fpass');

% Set the magUnits temporarily to 'linear' to get deviations
magUnits = get(d,'magUnits');
set(d,'magUnits','linear');
delta1 = get(d,'Dstop');
delta2 = get(d,'Dpass');
set(d,'magUnits',magUnits);

% [EOF]
