function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs.

%   Author(s): J. Schickler
%   Copyright 2006 The MathWorks, Inc.

specs.Fpass = this.Flow;
specs.Fstop = this.Fhigh;
specs.Apass = NaN;
specs.Astop = NaN;

% [EOF]
