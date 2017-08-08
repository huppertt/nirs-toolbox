function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fstop = (1-this.TransitionWidth)/2;
specs.Fpass = (1+this.TransitionWidth)/2;
specs.Astop = this.Astop;
specs.Apass = NaN;

% [EOF]
