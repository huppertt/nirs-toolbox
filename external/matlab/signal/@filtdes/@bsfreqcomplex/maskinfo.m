function cmd = maskinfo(hObj, d)
%MASKINFO Returns the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

Fn = getnyquist(d);

cmd{1}.frequency  = [-Fn get(d, 'Fpass1')];
cmd{1}.freqfcn    = 'wstop';
cmd{1}.filtertype = 'lowpass';
cmd{1}.tag        = 'pass1';

cmd{2}.frequency  = [get(d, 'Fstop1') get(d, 'Fstop2')];
cmd{2}.freqfcn    = 'wstop';
cmd{2}.filtertype = 'bandstop';

cmd{3}.frequency  = [get(d, 'Fpass2') get(d, 'Fpass3')];
cmd{3}.freqfcn    = 'wpass';
cmd{3}.filtertype = 'bandpass';

cmd{4}.frequency  = [get(d, 'Fstop3') get(d, 'Fstop4')];
cmd{4}.freqfcn    = 'wstop';
cmd{4}.filtertype = 'bandstop';

cmd{5}.frequency  = [get(d, 'Fpass4') getnyquist(d)];
cmd{5}.freqfcn    = 'wpass';
cmd{5}.filtertype = 'highpass';
cmd{5}.tag        = 'pass3';

% [EOF]
