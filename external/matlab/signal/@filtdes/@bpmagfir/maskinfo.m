function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

mu = get(d, 'magUnits');
if strcmpi(mu, 'db'),
    astop1 = get(d, 'Astop1');
    apass  = get(d, 'Apass');
    astop2 = get(d, 'Astop2');
    astopp = -max(astop1, astop2)-50;
else,
    astop1 = get(d, 'Dstop1');
    apass  = get(d, 'Dpass');
    astop2 = get(d, 'Dstop2');
    astopp = -max(astop1, astop2);
end

cmd{1}.magfcn     = 'stop';
cmd{1}.amplitude  = astop1;
cmd{1}.filtertype = 'highpass';
cmd{1}.magunits   = mu;
cmd{1}.tag        = 'stop1';

cmd{2}.magfcn     = 'pass';
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
