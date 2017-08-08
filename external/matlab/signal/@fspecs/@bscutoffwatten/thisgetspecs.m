function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs. This is used by FVTOOL for drawing the mask.

%   Copyright 2008 The MathWorks, Inc.

specs.Fstop1 = this.Fcutoff1;
specs.Fpass1 = this.Fcutoff1;
specs.Fpass2 = this.Fcutoff2;
specs.Fstop2 = this.Fcutoff2;
specs.Apass1 = this.Apass1;
specs.Astop = this.Astop;
specs.Apass2 = this.Apass2;

% [EOF]
