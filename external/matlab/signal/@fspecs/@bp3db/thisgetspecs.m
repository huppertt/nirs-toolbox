function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fstop1 = this.F3db1;
specs.Fpass1 = this.F3db1;
specs.Fpass2 = this.F3db2;
specs.Fstop2 = this.F3db2;
specs.Astop1 = NaN;
specs.Apass  = NaN;
specs.Astop2 = NaN;

% [EOF]
