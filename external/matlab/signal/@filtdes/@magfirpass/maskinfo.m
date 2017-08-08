function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

if isdb(d), apass = get(d, 'Apass');
else,       apass = get(d, 'Dpass'); end

cmd{1}.magfcn    = 'pass';
cmd{1}.amplitude = apass;
cmd{1}.magunits  = get(d, 'magUnits');

% [EOF]
