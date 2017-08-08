function cmd = maskinfo(hObj, d)
%MASKINFO Returns the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

Fn = getnyquist(d);

cmd{1}.frequency  = [-Fn get(d, 'Fstop1')];
cmd{1}.freqfcn    = 'wstop';
cmd{1}.filtertype = 'highpass';
cmd{1}.tag        = 'stop1';

cmd{2}.frequency  = [get(d, 'Fpass1') get(d, 'Fpass2')];
cmd{2}.freqfcn    = 'wpass';
cmd{2}.filtertype = 'bandpass';

cmd{3}.frequency  = [get(d, 'Fstop2') get(d, 'Fstop3')];
cmd{3}.freqfcn    = 'wstop';
cmd{3}.filtertype = 'bandstop';

cmd{4}.frequency  = [get(d, 'Fpass3') get(d, 'Fpass4')];
cmd{4}.freqfcn    = 'wpass';
cmd{4}.filtertype = 'bandpass';

cmd{5}.frequency  = [get(d, 'Fstop4') Fn];
cmd{5}.freqfcn    = 'wstop';
cmd{5}.filtertype = 'lowpass';
cmd{5}.tag        = 'stop2';

% [EOF]
