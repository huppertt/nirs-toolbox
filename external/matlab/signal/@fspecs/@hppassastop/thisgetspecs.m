function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fstop = this.Fpass;
specs.Fpass = this.Fpass;
specs.Astop = this.Astop;
specs.Apass = this.Apass;

% [EOF]
