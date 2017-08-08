function c = coeffs(this)
%COEFFS   Returns the coefficients.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% Add a field for each of the stages with that stage's coefficients.
for indx = 1:nstages(this)
    c.(sprintf('Stage%d', indx)) = coeffs(this.Stage(indx));
end

% [EOF]
