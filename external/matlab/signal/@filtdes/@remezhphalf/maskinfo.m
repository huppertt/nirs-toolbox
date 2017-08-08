function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

cmd = base_maskinfo(hObj, d);

cmd.bands{1}.freqfcn    = 'wpass';
cmd.bands{1}.frequency  = [get(d, 'Fpass') getnyquist(d)];
cmd.bands{1}.filtertype = 'highpass';

% [EOF]
