function minfo = measureinfo(this)
%MEASUREINFO   Return a structure of information for the measurements.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

minfo.Fstop = [];
minfo.Fpass = [];
minfo.Astop = this.Astop;
minfo.Apass = this.Apass;

% [EOF]
