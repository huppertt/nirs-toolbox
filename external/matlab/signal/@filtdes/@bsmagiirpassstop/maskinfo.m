function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

if isdb(d),
    apass = get(d, 'Apass');
    astop = get(d, 'Astop');
    astopp = -astop-50;
else,
    apass = get(d, 'Epass');
    astop = get(d, 'Estop');
    astopp = -astop;
end

mu = get(d, 'magUnits');
cmd{1}.magfcn     = 'cpass';
cmd{1}.amplitude  = apass;
cmd{1}.filtertype = 'lowpass';
cmd{1}.magunits   =  mu;
cmd{1}.astop      = astopp;
cmd{1}.tag        = 'pass1';

cmd{2}.magfcn     = 'stop';
cmd{2}.amplitude  = astop;
cmd{2}.filtertype = 'bandstop';
cmd{2}.magunits   =  mu;

cmd{3}.magfcn     = 'cpass';
cmd{3}.amplitude  = apass;
cmd{3}.filtertype = 'highpass';
cmd{3}.magunits   =  mu;
cmd{3}.astop      = astopp;
cmd{3}.tag        = 'pass2';

% [EOF]
