function minfo = measureinfo(this)
%MEASUREINFO   Return a structure of information for the measurements.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

astop = get(this, 'Astop');
apass = convertmagunits(convertmagunits(astop, 'db', 'linear', 'stop'), 'linear', 'db', 'pass');

minfo.Fpass = [];
minfo.Fstop = [];
minfo.Apass = apass;
minfo.Astop = astop;

% [EOF]
