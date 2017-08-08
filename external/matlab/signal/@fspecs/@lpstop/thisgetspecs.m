function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fpass = this.Fstop;
specs.Fstop = this.Fstop;
specs.Apass = NaN;
specs.Astop = this.Astop;

% [EOF]
