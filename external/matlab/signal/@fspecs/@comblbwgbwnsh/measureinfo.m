function minfo = measureinfo(this)
%MEASUREINFO   Return a structure of information for the measurements.

%   Copyright 2008 The MathWorks, Inc.

minfo.GBW    = this.GBW;
minfo.FilterOrder = this.NumPeaksOrNotches;
minfo.CombType = this.CombType;

% [EOF]
