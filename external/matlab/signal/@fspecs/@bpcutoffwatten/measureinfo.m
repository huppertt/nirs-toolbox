function minfo = measureinfo(this)
%MEASUREINFO   

%   Copyright 2008 The MathWorks, Inc.

minfo.Fstop1 = [];
minfo.Fcutoff1 = this.Fcutoff1;
minfo.Fpass1 = [];
minfo.Fpass2 = [];
minfo.Fcutoff2 = this.Fcutoff2;
minfo.Fstop2 = [];
minfo.Astop1 = this.Astop1;
minfo.Apass  = this.Apass;
minfo.Astop2 = this.Astop2;

% [EOF]
