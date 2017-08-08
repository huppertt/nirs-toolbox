function minfo = measureinfo(this)
%MEASUREINFO   Return a structure of information for the measurements.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

minfo.Fstop1 = [];
minfo.Fcutoff1 = [];
minfo.Fpass1 = this.Fpass1;
minfo.Fpass2 = this.Fpass2;
minfo.Fcutoff2 = [];
minfo.Fstop2 = [];
minfo.Astop1 = this.Astop1;
minfo.Apass  = this.Apass;
minfo.Astop2 = this.Astop2;

% [EOF]
