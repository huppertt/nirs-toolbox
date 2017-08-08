function lthoctave = set_bandsperoctave(this, lthoctave)
%SET_BANDSPEROCTAVE   PreSet function for the 'bandsperoctave' property.

%   Author(s): V. Pellissier
%   Copyright 2006 The MathWorks, Inc.

if isprop(this.CurrentSpecs, 'privBandsPerOctave')
    set(this.CurrentSpecs, 'privBandsPerOctave', lthoctave);
end


% [EOF]
