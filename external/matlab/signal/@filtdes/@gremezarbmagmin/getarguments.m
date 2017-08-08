function varargout = getarguments(h, d)
%GETARGUMENTS

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

F = get(d, 'FrequencyVector');
A = get(d, 'MagnitudeVector');
R = get(d, 'RippleVector');

if nargout == 1,
    varargout = {{F, A, R}};
else
    varargout = {F, A, R, {}};
end

% [EOF]
