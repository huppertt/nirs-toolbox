function specs = thisgetspecs(this)
%THISGETSPECS   

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

tw = (this.BWstop - (this.F3dB2 - this.F3dB1))/2;

specs.Fstop1 = this.F3dB1-tw;
specs.Fpass1 = this.F3dB1;
specs.Fpass2 = this.F3dB2;
specs.Fstop2 = this.F3dB2+tw;
specs.Astop1 = NaN;
specs.Apass  = NaN;
specs.Astop2 = NaN;

% [EOF]
