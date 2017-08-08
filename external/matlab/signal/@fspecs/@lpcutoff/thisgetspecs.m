function specs = thisgetspecs(this)
%THISGETSPECS   

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fpass = this.Fcutoff;
specs.Fstop = this.Fcutoff;
specs.Apass = NaN;
specs.Astop = NaN;

% [EOF]
