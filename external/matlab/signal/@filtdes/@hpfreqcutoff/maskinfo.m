function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

cmd = fc1_maskinfo(hObj, d);

cmd{1}.freqfcn = 'wstop';
cmd{2}.freqfcn = 'wpass';

cmd{1}.drawpatch    = false;
cmd{1}.drawfreqbars = false;

cmd{1}.filtertype = 'highpass';
cmd{2}.filtertype = 'highpass';

% [EOF]
