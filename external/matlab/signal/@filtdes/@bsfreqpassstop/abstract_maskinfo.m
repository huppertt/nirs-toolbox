function cmd = abstract_maskinfo(hObj, d)
%ABSTRACT_MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

cmd{1}.frequency = [0 get(d,'Fpass1')];
cmd{1}.freqfcn   = 'wpass';
cmd{1}.tag       = 'pass1';

cmd{2}.frequency = [get(d,'Fstop1'), get(d,'Fstop2')];
cmd{2}.freqfcn   = 'wstop';

cmd{3}.frequency = [get(d,'Fpass2'), getnyquist(d)];
cmd{3}.freqfcn   = 'wpass';
cmd{3}.tag       = 'pass2';

% [EOF]
