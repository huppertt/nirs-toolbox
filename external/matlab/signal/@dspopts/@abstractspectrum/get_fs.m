function Fs = get_fs(h, Fs) %#ok
%GETFS   Pre-Get Function for the Fs property.

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

if h.NormalizedFrequency,
    Fs = 'Normalized';
else
    Fs = h.privFs;
end

% [EOF]
