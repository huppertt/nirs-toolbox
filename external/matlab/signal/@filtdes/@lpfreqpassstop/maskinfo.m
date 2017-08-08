function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

cmd{1}.frequency  = [0 get(d,'Fpass')];
cmd{1}.filtertype = 'lowpass';
cmd{1}.freqfcn    = 'wpass';

cmd{2}.frequency  = [get(d,'Fstop') getnyquist(d)];
cmd{2}.filtertype = 'lowpass';
cmd{2}.freqfcn    = 'wstop';

% [EOF]
