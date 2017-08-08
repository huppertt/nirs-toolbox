function Hd = design(h, d)
%DESIGN Design the filter

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% Make sure we have the magnitude units in dB.
mu = get(d, 'MagUnits'); set(d, 'MagUnits', 'dB');
apass = get(d, 'Apass'); set(d, 'MagUnits', mu);

[b, a] = iirnotch(get(d, 'Fnotch'), getbandwidth(d), apass);
Hd     = dfilt.df2(b, a);

% [EOF]
