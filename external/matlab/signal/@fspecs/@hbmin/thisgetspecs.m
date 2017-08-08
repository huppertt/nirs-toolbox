function specs = thisgetspecs(this)
%THISGETSPECS   

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

specs.Fpass = (1-this.TransitionWidth)/2;
specs.Fstop = (1+this.TransitionWidth)/2;
specs.Apass = convertmagunits(convertmagunits(this.Astop, 'db', 'linear', 'stop'), ...
    'linear', 'db', 'pass');
specs.Astop = this.Astop;

% [EOF]
