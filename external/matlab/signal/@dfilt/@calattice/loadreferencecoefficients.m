function loadreferencecoefficients(this, s)
%LOADREFERENCECOEFFICIENTS   Load the reference coefficients.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

if s.version.number < 2
    set(this, 'Allpass1', s.Allpass1, ...
        'Allpass2', s.Allpass2, ...
        'Beta', s.Beta);
else
    set(this, 'Allpass1', s.refAllpass1, ...
        'Allpass2', s.refAllpass2, ...
        'Beta', s.refBeta);
end

% [EOF]
