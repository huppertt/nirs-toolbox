function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

tw = (this.F3dB1-this.F3dB2-this.BWstop)/2;

specs.Fpass1 = this.F3dB1;
specs.Fstop1 = this.F3dB1+tw;
specs.Fstop2 = this.F3dB2-tw;
specs.Fpass2 = this.F3dB2;
specs.Apass1 = NaN;
specs.Astop  = NaN;
specs.Apass2 = NaN;

% [EOF]
