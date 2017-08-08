function minfo = measureinfo(this)
%MEASUREINFO   Return a structure of information for the measurements.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

minfo.Fstop = this.Fstop;
minfo.Fpass = this.Fpass;
minfo.Astop = [];
minfo.Apass = [];

% [EOF]
