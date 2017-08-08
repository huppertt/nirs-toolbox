function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

if isdb(d),
    apass = get(d, 'Apass');
    astop = get(d, 'Astop2');
    astopp = -astop-50;
else
    apass = get(d, 'Dpass');
    astop = get(d, 'Dstop2');
    astopp = -astop;
end

cmd{1} = [];

cmd{2}.magfcn     = 'pass';
cmd{2}.amplitude  = apass;
cmd{2}.filtertype = 'bandpass';
cmd{2}.astop      = astopp;

cmd{3}.magfcn     = 'stop';
cmd{3}.amplitude  = astop;
cmd{3}.filtertype = 'lowpass';

% [EOF]
