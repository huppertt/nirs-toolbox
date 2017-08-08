function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

if isdb(d),
    apass = get(d, 'Apass1');
    astop = get(d, 'Astop');
    astopp = -astop-50;
else
    apass = get(d, 'Dpass1');
    astop = get(d, 'Dstop');
    astopp = -astop;
end

cmd{1}.magfcn     = 'pass';
cmd{1}.amplitude  = apass;
cmd{1}.filtertype = 'lowpass';
cmd{1}.astop      = astopp;

cmd{2}.magfcn     = 'stop';
cmd{2}.amplitude  = astop;
cmd{2}.filtertype = 'bandstop';

cmd{3} = [];

% [EOF]
