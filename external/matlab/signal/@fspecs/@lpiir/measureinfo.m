function minfo = measureinfo(this)
%MEASUREINFO   Return a structure of information for the measurements.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

minfo.Fpass = this.Fpass;
minfo.Fstop = this.Fstop;
minfo.Apass = [];
minfo.Astop = [];

% [EOF]
