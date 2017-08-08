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
cmd{2}.weight     = get(d, 'Wstop');

cmd{3}.magfcn     = 'wpass';
cmd{3}.filtertype = 'highpass';
cmd{3}.tag        = 'pass2';
cmd{3}.weight     = get(d, 'Wpass2');

% [EOF]
