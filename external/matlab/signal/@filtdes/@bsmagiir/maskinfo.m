function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

mu = get(d, 'magUnits');
if strcmpi(mu, 'db'),
    apass1 = get(d, 'Apass1');
    astop  = get(d, 'Astop');
    apass2 = get(d, 'Apass2');
    astopp = -astop-50;
else,
    apass1 = get(d, 'Epass1');
    astop  = get(d, 'Estop');
    apass2 = get(d, 'Epass2');
    astopp = -astop;
end

cmd{1}.magfcn     = 'cpass';
cmd{1}.amplitude  = apass1;
cmd{1}.filtertype = 'lowpass';
cmd{1}.magunits   = mu;
cmd{1}.tag        = 'pass1';
cmd{1}.astop      = astopp;

cmd{2}.magfcn     = 'stop';
cmd{2}.amplitude  = astop;
cmd{2}.filtertype = 'bandstop';
cmd{2}.magunits   = mu;

cmd{3}.magfcn     = 'cpass';
cmd{3}.amplitude  = apass2;
cmd{3}.filtertype = 'highpass';
cmd{3}.magunits   = mu;
cmd{3}.tag        = 'pass2';
cmd{3}.astop      = astopp;

% [EOF]
