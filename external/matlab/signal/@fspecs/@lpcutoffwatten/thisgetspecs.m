function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fpass = this.Fcutoff;
specs.Fstop = this.Fcutoff;
specs.Apass = this.Apass;
specs.Astop = this.Astop;

% [EOF]
