function minfo = measureinfo(this)
%MEASUREINFO   Return a structure of information for the measurements.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

minfo.Fpass1 = this.Fpass1;
minfo.Fcutoff1 = [];
minfo.Fstop1 = this.Fstop1;
minfo.Fstop2 = this.Fstop2;
minfo.Fcutoff2 = [];
minfo.Fpass2 = this.Fpass2;
minfo.Apass1 = this.Apass1;
minfo.Astop  = this.Astop;
minfo.Apass2 = this.Apass2;

% [EOF]
