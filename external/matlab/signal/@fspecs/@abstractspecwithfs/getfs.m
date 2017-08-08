function Fs = getfs(h,Fs)
%GETFS   Pre-Get Function for the Fs property.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

if h.NormalizedFrequency,
    Fs = 'Normalized';
else
    Fs = h.privFs;
end

% [EOF]
