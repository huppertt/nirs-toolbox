function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

if isdb(d),
    astop1 = get(d, 'Astop1');
    apass  = get(d, 'Apass');
    astop2 = get(d, 'Astop2');
    astopp = -max(astop1, astop2)-50;
else,
    astop1 = get(d, 'Estop1');
    apass  = get(d, 'Epass');
    astop2 = get(d, 'Estop2');
    astopp = -max(astop1, astop2);
end

mu = get(d, 'magUnits');

cmd{1}.magfcn     = 'stop';
cmd{1}.amplitude  = astop1;
cmd{1}.filtertype = 'highpass';
cmd{1}.magunits   = mu;
cmd{1}.tag        = 'stop1';

cmd{2}.magfcn     = 'cpass';
cmd{2}.amplitude  = apass;
cmd{2}.filtertype = 'bandpass';
cmd{2}.magunits   = mu;
cmd{2}.astop      = astopp;

cmd{3}.magfcn     = 'stop';
cmd{3}.amplitude  = astop2;
cmd{3}.filtertype = 'lowpass';
cmd{3}.magunits   = mu;
cmd{3}.tag        = 'stop2';


% [EOF]
