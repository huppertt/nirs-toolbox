function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

if isdb(d),
    apass = get(d, 'Apass');
    astop = max(60, apass+20);
else
    apass = get(d, 'Epass');
    astop = 0;
end

cmd{1}.magfcn     = 'cpass';
cmd{1}.amplitude  = apass;
cmd{1}.filtertype = 'lowpass';
cmd{1}.magunits   = get(d, 'magUnits');
cmd{1}.astop      = -astop;
cmd{1}.tag        = 'pass1';

cmd{2}            = [];

cmd{3}.magfcn     = 'cpass';
cmd{3}.amplitude  = apass;
cmd{3}.filtertype = 'highpass';
cmd{3}.magunits   = get(d, 'magUnits');
cmd{3}.astop      = -astop;
cmd{3}.tag        = 'pass2';

% [EOF]