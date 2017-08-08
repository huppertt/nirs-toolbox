function s = savereferencecoefficients(this)
%SAVEREFERENCECOEFFICIENTS   Save the reference coefficients.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

s.refsosMatrix   = get(this, 'refsosMatrix');
s.refScaleValues = get(this, 'refScaleValues');

% [EOF]
