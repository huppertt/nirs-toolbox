function delay = set_delay(this, delay)
%SET_DELAY   PreSet function for the 'delay' property.

%   Author(s): V. Pellissier
%   Copyright 2005-2006 The MathWorks, Inc.

Fs = this.Fs;
if delay<0,
    error(message('signal:fdesign:fracdelay:set_delay:InvalidFracDelay'));
end

if isprop(this.CurrentSpecs, 'privFracDelay')
    set(this.CurrentSpecs, 'privFracDelay', delay);
end

% [EOF]
