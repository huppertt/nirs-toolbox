function s = savereferencecoefficients(this)
%SAVEREFERENCECOEFFICIENTS   Save the reference coefficients.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

s.refAllpass1 = get(this, 'refAllpass1');
s.refAllpass2 = get(this, 'refAllpass2');
s.refBeta     = get(this, 'refBeta');

% [EOF]
