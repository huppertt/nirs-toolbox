function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

mu = get(d, 'magUnits');
if isdb(d), astop = get(d, 'Astop');
else,       astop = get(d, 'Estop'); end

cmd{1}.magfcn     = 'stop';
cmd{1}.amplitude  = astop;
cmd{1}.filtertype = 'highpass';
cmd{1}.magunits   = mu;
cmd{1}.tag        = 'stop1';
cmd{1}.drawpatch  = true;

cmd{2} = [];

cmd{3}.magfcn     = 'stop';
cmd{3}.amplitude  = astop;
cmd{3}.filtertype = 'lowpass';
cmd{3}.magunits   = mu;
cmd{3}.tag        = 'stop2';
cmd{3}.drawpatch  = true;

% [EOF]
