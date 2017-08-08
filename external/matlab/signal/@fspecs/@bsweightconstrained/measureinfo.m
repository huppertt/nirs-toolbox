function minfo = measureinfo(this)
%MEASUREINFO Return a structure of information for the measurements.

%   Copyright 2011 The MathWorks, Inc.

minfo.Fpass1 = this.Fpass1;
minfo.Fcutoff1 = [];
minfo.Fstop1 = this.Fstop1;
minfo.Fstop2 = this.Fstop2;
minfo.Fcutoff2 = [];
minfo.Fpass2 = this.Fpass2;
if this.Passband1Constrained
  minfo.Apass1 = this.Apass1;
else
  minfo.Apass1 = [];
end
if this.StopbandConstrained
  minfo.Astop = this.Astop;
else
  minfo.Astop = [];
end
if this.Passband2Constrained
  minfo.Apass2  = this.Apass2;
else  
  minfo.Apass2  = [];
end

% [EOF]
