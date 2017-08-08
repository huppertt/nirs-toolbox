function band = set_band(this, band)
%SET_BAND   PreSet function for the 'band' property.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

if band <= 1,
    error(message('signal:fdesign:abstracttypewband:set_band:invalidBand'));
end

set(this, 'privBand', band);

% Make sure that the current specifications are up to date based on the new
% value of the band property.
updatecurrentspecs(this);

% If the band property exists on the new specifications, set it.
if isprop(this.CurrentSpecs, 'Band')
    set(this.CurrentSpecs, 'Band', band);
end

% [EOF]
