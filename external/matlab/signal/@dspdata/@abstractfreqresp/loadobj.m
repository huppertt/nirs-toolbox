function this = loadobj(s)
%LOADOBJ   Load this object.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

this = feval(s.class);

normalizefreq(this, s.NormalizedFrequency);

set(this, ...
    'privFs',      s.Fs, ...
    'Data',        s.Data, ...
    'Frequencies', s.Frequencies, ...
    'CenterDC',    s.CenterDC, ...
    'Metadata',    s.Metadata);

thisloadobj(this, s);

% [EOF]
