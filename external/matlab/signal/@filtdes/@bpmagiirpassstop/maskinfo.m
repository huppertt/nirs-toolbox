function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

mu = get(d, 'magUnits');
if isdb(d),
    apass = get(d, 'Apass');
    astop = get(d, 'Astop');
    astopp = -astop-50;
else,
    apass = get(d, 'Epass');
    astop = get(d, 'Estop');
    astopp = -astop;
end

cmd{1}.magfcn     = 'stop';
cmd{1}.amplitude  = astop;
cmd{1}.filtertype = 'highpass';
cmd{1}.magunits   = mu;
cmd{1}.tag        = 'stop1';

cmd{2}.magfcn     = 'cpass';
cmd{2}.amplitude  = apass;
cmd{2}.filtertype = 'bandpass';
cmd{2}.magunits   = mu;
cmd{2}.astop      = astopp;

cmd{3}.magfcn     = 'stop';
cmd{3}.amplitude  = astop;
cmd{3}.filtertype = 'lowpass';
cmd{3}.magunits   = mu;
cmd{3}.tag        = 'stop2';

% [EOF]
