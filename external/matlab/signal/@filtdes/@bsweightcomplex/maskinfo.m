function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

cmd{1}.magfcn     = 'wpass';
cmd{1}.filtertype = 'lowpass';
cmd{1}.tag        = 'pass1';
cmd{1}.weight     = get(d, 'Wpass1');

cmd{2}.magfcn     = 'wstop';
cmd{2}.filtertype = 'bandstop';
cmd{2}.weight     = get(d, 'Wstop1');

cmd{3}.magfcn     = 'wpass';
cmd{3}.filtertype = 'bandpass';
cmd{3}.weight     = get(d, 'Wpass2');

cmd{4}.magfcn     = 'wstop';
cmd{4}.filtertype = 'bandstop';
cmd{4}.weight     = get(d, 'Wstop2');

cmd{5}.magfcn     = 'wpass';
cmd{5}.filtertype = 'highpass';
cmd{5}.tag        = 'pass2';
cmd{5}.weight     = get(d, 'Wpass3');

% [EOF]
