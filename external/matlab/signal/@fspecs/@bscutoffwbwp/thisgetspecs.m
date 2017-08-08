function specs = thisgetspecs(this)
%THISGETSPECS   

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

tw = (this.BWpass - (this.F3dB2-this.F3dB1))/2;

specs.Fpass1 = this.F3dB1-tw;
specs.Fstop1 = this.F3dB1;
specs.Fstop2 = this.F3dB2;
specs.Fpass2 = this.F3dB2+tw;
specs.Apass1 = NaN;
specs.Astop  = NaN;
specs.Apass2 = NaN;

% [EOF]
