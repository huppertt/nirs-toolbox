function minfo = measureinfo(this)
%MEASUREINFO Return a structure of information for the measurements.

%   Copyright 2011 The MathWorks, Inc.

minfo.Fpass = this.Fpass;
minfo.Fstop = this.Fstop;
minfo.Apass = this.Apass;
minfo.Astop = [];

% [EOF]
