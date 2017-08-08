function specs = thisgetspecs(this)
%THISGETSPECS   

%   Copyright 2011 The MathWorks, Inc.

specs.Fpass1 = this.Fpass1;
specs.Fstop1 = this.Fstop1;
specs.Fstop2 = this.Fstop2;
specs.Fpass2 = this.Fpass2;
if this.Passband1Constrained
    specs.Apass1  = this.Apass1;
else
  specs.Apass1  = NaN;
end
if this.StopbandConstrained
  specs.Astop = this.Astop;
else
  specs.Astop = NaN;
end
if this.Passband2Constrained
    specs.Apass2  = this.Apass2;
else
  specs.Apass2  = NaN;
end

% [EOF]
