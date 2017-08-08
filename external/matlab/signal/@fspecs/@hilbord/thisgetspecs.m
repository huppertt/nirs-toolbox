function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

TWn = this.TransitionWidth/2;
specs.Fstop = 1-TWn;
specs.Fpass = TWn;
specs.Apass = nan;
specs.Astop = nan;

% [EOF]
