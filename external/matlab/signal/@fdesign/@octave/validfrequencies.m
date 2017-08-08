function f = validfrequencies(this)
%VALIDFREQUENCIES  Return the valid values for the 'CenterFreq' property.

%   Author(s): V. Pellissier
%   Copyright 2006 The MathWorks, Inc.

f = getvalidcenterfrequencies(this.Currentspecs);
if this.NormalizedFrequency,
    f = f/24000;
end

% [EOF]
