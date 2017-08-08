function specs = thisgetspecs(this)
%THISGETSPECS   Get specifications.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

specs.Fpass  = this.Fpass;
specs.Fstop = nan;
specs.Apass = nan;
specs.Astop = this.Astop;

% [EOF]
