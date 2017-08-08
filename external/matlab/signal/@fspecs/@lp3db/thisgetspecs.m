function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fpass = this.F3dB;
specs.Fstop = this.F3dB;
specs.Apass = NaN;
specs.Astop = NaN;

% [EOF]
