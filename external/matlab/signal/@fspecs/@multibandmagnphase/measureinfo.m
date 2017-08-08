function minfo = measureinfo(this)
%MEASUREINFO   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

[F, H] = getmask(this);
minfo.Frequencies = F;
minfo.FreqResponse = H;


% [EOF]
