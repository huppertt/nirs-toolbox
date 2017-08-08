function minfo = measureinfo(this)
%MEASUREINFO Return a structure of information for the measurements.

%   Copyright 2011 The MathWorks, Inc.

minfo.Fpass = [];
minfo.Fcutoff = [];
minfo.Fstop = this.Fstop;
minfo.Apass = this.Apass;
minfo.Astop = this.Astop;

minfo.FrequencyFactor = this.FrequencyFactor;
minfo.Power           = this.Power;

% [EOF]
