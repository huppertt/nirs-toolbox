function cmd = fs1_maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

cmd{1}.frequency = [0 get(d, 'Fstop')];
cmd{2}.frequency = [get(d, 'Fstop') getnyquist(d)];

% [EOF]
