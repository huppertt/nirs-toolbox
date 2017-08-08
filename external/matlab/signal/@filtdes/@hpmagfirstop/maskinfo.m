function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

if isdb(d),
    astop = get(d, 'Astop');
else
    astop = get(d, 'Dstop');
end

cmd{1}.magfcn     = 'stop';
cmd{1}.amplitude  = astop;
cmd{1}.filtertype = 'highpass';

cmd{2} = {};

% [EOF]
