function [Dpass, Dstop] = getdesignspecs(h, d)
%GETDESIGNSPECS Returns the specs for the design

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Get mag specs
magUnits = get(d,'magUnits');
set(d,'magUnits','linear'); % Set temporarily to linear
Dpass = get(d,'Dpass');
Dstop = get(d,'Dstop');
set(d,'magUnits',magUnits); % Reset units

% [EOF]
