function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fpass1 = this.F3db1;
specs.Fstop1 = this.F3db1;
specs.Fstop2 = this.F3db2;
specs.Fpass2 = this.F3db2;
specs.Apass1 = NaN;
specs.Astop  = NaN;
specs.Apass2 = NaN;

% [EOF]
