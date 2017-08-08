function wf = whichframes(h)
%WHICHFRAMES

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

wf = super_whichframes(h);

% Replace the textOptions frame with the cremezoptsframe
indx = strcmpi('siggui.textOptionsFrame', {wf.constructor});

if ismethod(h.ResponseTypeSpec, 'getoptsframe'),
    wf(indx) = getoptsframe(h.ResponseTypeSpecs);
else
    wf(indx).constructor = 'siggui.cremezoptsframe';
    wf(indx).setops      = {};
end

% [EOF]
