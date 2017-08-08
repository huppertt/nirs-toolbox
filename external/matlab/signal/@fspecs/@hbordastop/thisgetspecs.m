function specs = thisgetspecs(this)
%THISGETSPECS   Get the specs.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fpass = 1/2;
specs.Fstop = 1/2;
specs.Apass = convertmagunits(convertmagunits(this.Astop, 'db', 'linear', 'stop'), ...
    'linear', 'db', 'pass');
specs.Astop = this.Astop;

% [EOF]
