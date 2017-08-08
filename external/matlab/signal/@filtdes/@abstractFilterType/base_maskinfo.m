function cmd = base_maskinfo(hObj, d)
%BASE_MASKINFO Return the generic information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

cmd.fs        = getnyquist(d)*2;
cmd.frequnits = get(d, 'freqUnits');
cmd.response  = maskresponse(hObj);

if isprop(d, 'magUnits'),
    cmd.magunits = get(d, 'magUnits');
else
    cmd.magunits = 'weights';
end

% [EOF]
