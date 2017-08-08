function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fstop = this.F3db;
specs.Fpass = this.F3db;
specs.Astop = NaN;
specs.Apass = this.Apass;

% [EOF]
