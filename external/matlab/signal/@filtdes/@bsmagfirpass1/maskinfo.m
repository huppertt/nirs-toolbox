function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

if isdb(d),
    apass = get(d, 'Apass1');
else
    apass = get(d, 'Dpass1');
end

cmd{1}.magfcn     = 'pass';
cmd{1}.amplitude  = apass;
cmd{1}.filtertype = 'lowpass';

cmd{2} = [];
cmd{3} = [];

% [EOF]
