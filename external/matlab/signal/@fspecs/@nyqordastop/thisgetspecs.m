function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fpass = 1/this.Band;
specs.Fstop = 1/this.Band;
specs.Apass = NaN;
specs.Astop = this.Astop;

% [EOF]
