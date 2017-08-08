function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

if isdb(d), astop = get(d, 'Astop');
else,       astop = get(d, 'Estop'); end

cmd{1}.magfcn     = 'stop';
cmd{1}.amplitude  = astop;
cmd{1}.filtertype = 'highpass';
cmd{1}.magunits   = get(d, 'magUnits');
cmd{1}.drawpatch  = true;

cmd{2}            = [];

% [EOF]
