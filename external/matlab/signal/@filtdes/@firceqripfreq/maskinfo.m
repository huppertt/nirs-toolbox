function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

propname = determine_dynamicprop(d,...
    get(d, 'freqSpecType'),set(d, 'freqSpecType'));
Fc = get(d, propname);

cmd{1}.frequency  = [0 Fc];
cmd{2}.frequency  = [Fc getnyquist(d)];

% [EOF]
