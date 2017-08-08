function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fstop = this.Fcutoff;
specs.Fpass = this.Fcutoff;
specs.Astop = NaN;
specs.Apass = NaN;

% [EOF]
