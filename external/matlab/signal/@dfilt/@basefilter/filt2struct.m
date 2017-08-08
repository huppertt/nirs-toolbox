function s = filt2struct(this)
%FILT2STRUCT   Return a structure representation of the object.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% Add a field for the class.
s.class = class(this);

% Add the coefficients.
s       = setstructfields(s, coeffs(this));

% [EOF]
