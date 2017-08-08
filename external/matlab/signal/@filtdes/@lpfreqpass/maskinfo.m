function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

cmd = fp1_maskinfo(hObj, d);

cmd{1}.freqfcn = 'wpass';
cmd{2}.freqfcn = 'wstop';

cmd{2}.drawpatch    = false;
cmd{2}.drawfreqbars = false;

cmd{1}.filtertype = 'lowpass';
cmd{2}.filtertype = 'lowpass';

% [EOF]
