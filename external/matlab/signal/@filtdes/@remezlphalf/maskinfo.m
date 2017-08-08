function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

cmd = base_maskinfo(hObj, d);

cmd.bands{1}.magfcn    = 'wpass';
cmd.bands{1}.frequency  = [0 get(d, 'Fpass')];
cmd.bands{1}.filtertype = 'lowpass';

% [EOF]
