function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

if isdb(d),
    astop = get(d, 'Astop1');
    apass = get(d, 'Apass');
    astopp = -astop-50;
else
    astop = get(d, 'Dstop1');
    apass = get(d, 'Dpass');
    astopp = -astop;
end

cmd{1}.magfcn     = 'stop';
cmd{1}.amplitude  = astop;
cmd{1}.filtertype = 'highpass';

cmd{2}.magfcn     = 'pass';
cmd{2}.amplitude  = apass;
cmd{2}.filtertype = 'bandpass';
cmd{2}.astop      = astopp;

cmd{3} = [];

% [EOF]
