function minfo = measureinfo(this)
%MEASUREINFO   

%   Copyright 2011 The MathWorks, Inc.

[F, H] = getmask(this);
minfo.Frequencies = F;
minfo.FreqResponse = H;


% [EOF]
