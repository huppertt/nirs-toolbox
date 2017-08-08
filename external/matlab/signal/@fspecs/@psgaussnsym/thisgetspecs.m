function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs.

%   Copyright 2008 The MathWorks, Inc.

specs.Fpass = 2*this.BT/this.SamplesPerSymbol;
specs.Fstop = NaN;
specs.Apass = NaN;
specs.Astop = NaN;

% [EOF]
