function y = nonzerosign(x)
%

%NONZEROSIGN Signum function excluding zero
%   For each element of X, NONZEROSIGN(X) returns 1 if the element
%   is greater than or equal to zero and -1 if it is
%   less than zero. NONZEROSIGN differs from SIGN in that NONZEROSIGN(0)
%   returns 1, while SIGN(0) returns 0.
%

%   Copyright 2008 The MathWorks, Inc.

y = ones(size(x));

y(x < 0) = -1;