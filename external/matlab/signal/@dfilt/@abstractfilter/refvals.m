function v = refvals(this)
%REFVALS   Return the reference values.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% Default refvals just gets the coefficients
c = coefficientnames(this);
for indx = 1:length(c)
    v{indx} = get(this, c{indx});
end

% [EOF]
