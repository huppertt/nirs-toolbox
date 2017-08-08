function coeffs = get_coeffs(this, coeffs)
%GET_COEFFS   PreGet function for the 'coeffs' property.

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

% Use getnumerator since it simply casts to double
coeffs = getnumerator(this.filterquantizer, this.privallpasscoeffs);

% [EOF]
