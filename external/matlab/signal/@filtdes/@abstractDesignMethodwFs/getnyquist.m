function fn = getnyquist(d)
%GETNYQUIST Returns the nyquist frequency.

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

if strncmpi(get(d, 'freqUnits'), 'normalized', 10),
    fn = 1;
else
    fn = get(d, 'Fs')/2;
end

% [EOF]
