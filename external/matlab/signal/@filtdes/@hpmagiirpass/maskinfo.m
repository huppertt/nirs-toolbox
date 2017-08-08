function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

if isdb(d), apass = get(d, 'Apass');
else,       apass = get(d, 'Epass'); end

cmd{1}            = [];

cmd{2}.magfcn     = 'cpass';
cmd{2}.amplitude  = apass;
cmd{2}.filtertype = 'highpass';
cmd{2}.magunits   = get(d, 'magUnits');

% [EOF]
