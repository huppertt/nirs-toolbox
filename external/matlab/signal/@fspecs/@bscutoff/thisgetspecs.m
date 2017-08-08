function specs = thisgetspecs(this)
%THISGETSPECS   

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fpass1 = this.Fcutoff1;
specs.Fstop1 = this.Fcutoff1;
specs.Fstop2 = this.Fcutoff2;
specs.Fpass2 = this.Fcutoff2;
specs.Apass1 = NaN;
specs.Astop  = NaN;
specs.Apass2 = NaN;

% [EOF]
