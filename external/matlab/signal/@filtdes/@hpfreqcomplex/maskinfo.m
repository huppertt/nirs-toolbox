function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

cmd = abstract_maskinfo(hObj, d);

cmd{1}.frequency(1) = -getnyquist(d);

% [EOF]
