function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs.

%   Copyright 2008 The MathWorks, Inc.

specs.Fpass = (1-this.RolloffFactor)/this.SamplesPerSymbol;
specs.Fstop = (1+this.RolloffFactor)/this.SamplesPerSymbol;
specs.Apass = NaN;
specs.Astop = this.Astop;

% [EOF]
