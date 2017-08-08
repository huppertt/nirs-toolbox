function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

if isdb(d),
    astop = get(d, 'Astop');
    apass = get(d, 'Apass2');
    astopp = -astop-50;
else
    astop = get(d, 'Dstop');
    apass = get(d, 'Dpass2');
    astopp = -astop;
end

cmd{1} = [];

cmd{2}.magfcn     = 'stop';
cmd{2}.amplitude  = astop;
cmd{2}.filtertype = 'bandstop';

cmd{3}.magfcn     = 'pass';
cmd{3}.amplitude  = apass;
cmd{3}.filtertype = 'highpass';
cmd{3}.astop      = astopp;

% [EOF]
