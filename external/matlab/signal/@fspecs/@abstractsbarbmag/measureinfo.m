function minfo = measureinfo(this)
%MEASUREINFO Return a structure of information for the measurements.

%   Copyright 2005-2011 The MathWorks, Inc.

[F, A] = getmask(this);
minfo.Frequencies = F;
minfo.Amplitudes = A;

% [EOF]
