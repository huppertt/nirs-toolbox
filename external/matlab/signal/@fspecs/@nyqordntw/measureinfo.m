function minfo = measureinfo(this)
%MEASUREINFO   Return a structure of information for the measurements.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

tw = get(this, 'TransitionWidth');

if this.NormalizedFrequency, Fs = 2;
else,                        Fs = this.Fs; end

band = get(this, 'Band');

minfo.Fpass = Fs/2/band-tw/2;
minfo.Fstop = Fs/2/band+tw/2;
minfo.Apass = [];
minfo.Astop = [];

% [EOF]
