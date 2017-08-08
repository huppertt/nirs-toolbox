function cmd = abstract_maskinfo(hObj, d)
%MASKINFO Returns the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

cmd{1}.frequency  = [0 get(d, 'Fstop1')];
cmd{1}.freqfcn    = 'wstop';
cmd{1}.filtertype = 'highpass';
cmd{1}.tag        = 'stop1';

cmd{2}.frequency  = [get(d, 'Fpass1') get(d, 'Fpass2')];
cmd{2}.freqfcn    = 'wpass';
cmd{2}.filtertype = 'bandpass';

cmd{3}.frequency  = [get(d, 'Fstop2') getnyquist(d)];
cmd{3}.freqfcn    = 'wstop';
cmd{3}.filtertype = 'lowpass';
cmd{3}.tag        = 'stop2';

% [EOF]
