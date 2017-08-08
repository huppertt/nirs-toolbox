function cmd = maskinfo(h,d)
%MASKINFO

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

cmd = firceqriphp_maskinfo(h,d);

cmd.bands{2}.magfcn     = 'invsinc';
cmd.bands{2}.frequency  = cmd.bands{2}.frequency(1);
cmd.bands{2}.freqfactor = get(d, 'InvSincFreqFactor');
cmd.bands{2}.sincpower  = get(d, 'InvSincPower');

% [EOF]
