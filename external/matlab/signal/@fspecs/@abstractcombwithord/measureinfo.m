function minfo = measureinfo(this)
%MEASUREINFO   Return a structure of information for the measurements.

%   Copyright 2008 The MathWorks, Inc.

minfo.GBW    = 10*log10(.5);
minfo.FilterOrder = this.FilterOrder;
minfo.CombType = this.CombType;

% [EOF]
