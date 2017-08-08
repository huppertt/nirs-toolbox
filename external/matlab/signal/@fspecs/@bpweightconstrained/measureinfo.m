function minfo = measureinfo(this)
%MEASUREINFO Return a structure of information for the measurements.

%   Copyright 2011 The MathWorks, Inc.

minfo.Fstop1 = this.Fstop1;
minfo.Fcutoff1 = [];
minfo.Fpass1 = this.Fpass1;
minfo.Fpass2 = this.Fpass2;
minfo.Fcutoff2 = [];
minfo.Fstop2 = this.Fstop2;
if this.Stopband1Constrained
  minfo.Astop1 = this.Astop1;
else
  minfo.Astop1 = [];
end
if this.PassbandConstrained
  minfo.Apass  = this.Apass;
else  
  minfo.Apass  = [];
end
if this.Stopband2Constrained
  minfo.Astop2 = this.Astop2;
else
  minfo.Astop2 = [];
end

% [EOF]
