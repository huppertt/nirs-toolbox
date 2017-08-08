function specs = thisgetspecs(this)
%THISGETSPECS   

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fpass1 = this.Fpass1;
specs.Fstop1 = this.Fpass1;
specs.Fstop2 = this.Fpass2;
specs.Fpass2 = this.Fpass2;
specs.Apass1 = this.Apass;
specs.Astop  = NaN;
specs.Apass2 = this.Apass;

% [EOF]
