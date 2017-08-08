function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fpass = this.Fpass;
specs.Fstop = this.Fpass;
specs.Apass = this.Apass;
specs.Astop = NaN;

% [EOF]
