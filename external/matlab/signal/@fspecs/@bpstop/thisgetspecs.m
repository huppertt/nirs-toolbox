function specs = thisgetspecs(this)
%THISGETSPECS   

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fstop1 = this.Fstop1;
specs.Fpass1 = this.Fstop1;
specs.Fpass2 = this.Fstop2;
specs.Fstop2 = this.Fstop2;
specs.Astop1 = this.Astop;
specs.Apass  = NaN;
specs.Astop2 = this.Astop;

% [EOF]
