function specs = thisgetspecs(this)
%THISGETSPECS   

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fstop1 = this.Fpass1;
specs.Fpass1 = this.Fpass1;
specs.Fpass2 = this.Fpass2;
specs.Fstop2 = this.Fpass2;
specs.Astop1 = this.Astop1;
specs.Apass  = this.Apass;
specs.Astop2 = this.Astop2;

% [EOF]
