function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

if isdb(d),
    apass = get(d, 'Apass2');
else
    apass = get(d, 'Dpass2');
end

cmd{1} = [];
cmd{2} = [];

cmd{3}.magfcn     = 'pass';
cmd{3}.amplitude  = apass;
cmd{3}.filtertype = 'highpass';

% [EOF]
