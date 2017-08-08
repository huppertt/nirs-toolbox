function [Fpass, Dpass] = getdesignspecs(h, d)
%GETDESIGNSPECS

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% get passband edge, it has been prenormalized
Fpass = get(d,'Fpass');

% Set the magUnits temporarily to 'linear' to get deviations
magUnits = get(d,'magUnits');
set(d,'magUnits','linear');
Dpass = get(d,'Dpass');
set(d,'magUnits',magUnits);

% [EOF]
