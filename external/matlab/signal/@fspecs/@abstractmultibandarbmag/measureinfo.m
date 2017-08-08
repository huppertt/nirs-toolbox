function minfo = measureinfo(this)
%MEASUREINFO   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

[F, A] = getmask(this);
minfo.Frequencies = F;
minfo.Amplitudes = A;


% [EOF]
