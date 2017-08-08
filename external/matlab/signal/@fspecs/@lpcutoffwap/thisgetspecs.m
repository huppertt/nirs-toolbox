function specs = thisgetspecs(this)
%THISGETSPECS   

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fpass = this.F3dB;
specs.Fstop = this.F3dB;
specs.Apass = this.Apass;
specs.Astop = NaN;

% [EOF]
