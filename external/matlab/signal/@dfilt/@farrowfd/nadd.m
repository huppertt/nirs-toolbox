function n = nadd(this)
%NADD Returns the number of adders  

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

c = coefficients(this);
coeffs = c{1};
nphases = size(coeffs,2);

n = nmult(this,false,false)-nphases;
if n<0, n = 0; end

% [EOF]
