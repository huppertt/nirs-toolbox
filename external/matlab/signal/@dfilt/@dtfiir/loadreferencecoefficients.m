function loadreferencecoefficients(this, s)
%LOADREFERENCECOEFFICIENTS   Load the reference coefficients.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

if s.version.number < 2
    set(this, 'Numerator', s.Numerator, ...
        'Denominator', s.Denominator);
else
    set(this, 'Numerator', s.refnum, ...
        'Denominator', s.refden);
end

% [EOF]
