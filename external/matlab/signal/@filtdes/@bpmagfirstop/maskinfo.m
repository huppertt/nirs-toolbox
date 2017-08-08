function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

if isdb(d),
    astop1 = get(d, 'Astop1');
    astop2 = get(d, 'Astop2');
else
    astop1 = get(d, 'Dstop1');
    astop2 = get(d, 'Dstop2');
end

cmd{1}.magfcn     = 'stop';
cmd{1}.amplitude  = astop1;
cmd{1}.filtertype = 'highpass';

cmd{2} = [];

cmd{3}.magfcn     = 'stop';
cmd{3}.amplitude  = astop2;
cmd{3}.filtertype = 'lowpass';

% [EOF]
