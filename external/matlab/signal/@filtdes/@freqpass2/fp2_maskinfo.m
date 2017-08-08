function cmd = fp2_maskinfo(hObj, d)
%FP2_MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

cmd{1}.frequency  = [0 get(d, 'Fpass1')];
cmd{2}.frequency  = [get(d, 'Fpass1') get(d, 'Fpass2')];
cmd{3}.frequency  = [get(d, 'Fpass2') getnyquist(d)];

% [EOF]
