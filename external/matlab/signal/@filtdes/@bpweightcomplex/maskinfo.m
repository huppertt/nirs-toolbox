function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

cmd{1}.magfcn     = 'wstop';
cmd{1}.filtertype = 'highpass';
cmd{1}.tag        = 'stop1';
cmd{1}.weight     = get(d, 'Wstop1');

cmd{2}.magfcn     = 'wpass';
cmd{2}.filtertype = 'bandpass';
cmd{2}.weight     = get(d, 'Wpass1');

cmd{3}.magfcn     = 'wstop';
cmd{3}.filtertype = 'bandstop';
cmd{3}.weight     = get(d, 'wstop2');

cmd{4}.magfcn     = 'wpass';
cmd{4}.filtertype = 'bandpass';
cmd{4}.weight     = get(d, 'Wpass2');

cmd{5}.magfcn     = 'wstop';
cmd{5}.filtertype = 'lowpass';
cmd{5}.tag        = 'stop2';
cmd{5}.weight     = get(d, 'Wstop3');

% [EOF]
