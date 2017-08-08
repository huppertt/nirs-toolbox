function s = savereferencecoefficients(this)
%SAVEREFERENCECOEFFICIENTS   Save the reference coefficients.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

s.refA = get(this, 'refA');
s.refB = get(this, 'refB');
s.refC = get(this, 'refC');
s.refD = get(this, 'refD');

% [EOF]
