function [Fstop, Fpass, Astop, Apass] = getdesignspecs(~, d)
%GETDESIGNSPECS Returns the specs for the design

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.

% Get frequency specs, they have been prenormalized
Fstop = get(d,'Fstop');
Fpass = get(d,'Fpass');

% Set the magUnits temporarily to 'dB' to get attenuations
magUnits = get(d,'magUnits');
%do some error checking on values of estop and epass
if strcmpi(magUnits, 'Squared')
    %Now check if Estop and Epass are greater than zero but less than 1.
    %If not, throw an error
    if d.Estop <= 0 || d.Estop >= 1 || d.Epass <= 0 || d.Epass >= 1
        error(message('signal:fdatool:EpassEstopLessThanZeroError'));
    end
    %Now check that Estop > Epass as the pass area should be higher than
    %the stop limit
    if d.Estop > d.Epass
        error(message('signal:fdatool:EpassGreaterThanEstopError'));
    end
end

% Set the magUnits temporarily to 'dB' to get attenuations
magUnits = get(d,'magUnits');
set(d,'magUnits','dB');
Astop = get(d,'Astop');
Apass = get(d,'Apass');

%Now check if Astop and Apass are greater than zero
%If not, throw an error
if d.Astop <= 0 || d.Apass <= 0
    error(message('signal:fdatool:ApassAstopLessThanZeroError'));
end
%Now check that Apass > Astop as the pass band should be smaller than
%the stop band
if d.Astop < d.Apass
    error(message('signal:fdatool:AstopLessThanApassError'));
end

set(d,'magUnits',magUnits);

% [EOF]
