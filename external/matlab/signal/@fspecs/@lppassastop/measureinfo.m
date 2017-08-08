function minfo = measureinfo(this)
%MEASUREINFO   Return a structure of information for the measurements.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

minfo.Fpass = this.Fpass;
minfo.Fstop = [];
minfo.Apass = this.Apass;
minfo.Astop = this.Astop;

% [EOF]
