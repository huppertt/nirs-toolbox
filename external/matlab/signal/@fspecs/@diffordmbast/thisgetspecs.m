function specs = thisgetspecs(this)
%THISGETSPECS Get the specs.

%   Copyright 2011 The MathWorks, Inc.

specs.Fpass = this.Fpass;
specs.Fstop = this.Fstop;
specs.Apass = nan;
specs.Astop = this.Astop;

% [EOF]
