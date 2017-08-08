function wf = whichframes(h)
%WHICHFRAMES

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

wf = anom_whichframes(h);

indx = find(strcmpi({wf.constructor}, 'siggui.textOptionsFrame'));

if isempty(indx), indx = length(wf)+1; end

wf(indx).constructor = 'siggui.ifiroptsframe';
wf(indx).setops      = {};

% [EOF]
