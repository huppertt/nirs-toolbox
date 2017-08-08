function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

if isdb(d),
    apass = get(d, 'Apass');
    astop = get(d, 'Astop');
    astopp = -astop-50;
else,
    apass = get(d, 'Dpass');
    astop = get(d, 'Dstop');
    astopp = -astop;
end

cmd{1}.magfcn     = 'stop';
cmd{1}.amplitude  = astop;
cmd{1}.filtertype = 'highpass';
cmd{1}.magunits   = get(d, 'magUnits');

cmd{2}.magfcn     = 'pass';
cmd{2}.amplitude  = apass;
cmd{2}.filtertype = 'highpass';
cmd{2}.magunits   = get(d, 'magUnits');
cmd{2}.astop      = astopp;

% [EOF]
