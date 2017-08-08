function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

if isdb(d),
    astop = get(d, 'Astop');
else
    astop = get(d, 'Dstop');
end

cmd{1} = [];

cmd{2}.magfcn     = 'stop';
cmd{2}.amplitude  = astop;
cmd{2}.filtertype = 'bandstop';

cmd{3} = [];

% [EOF]
