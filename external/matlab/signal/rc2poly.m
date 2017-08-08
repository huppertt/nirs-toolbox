function [a,efinal]=rc2poly(kr,R0)
%RC2POLY Convert reflection coefficients to prediction polynomial (step-up).
%   A = RC2POLY(K) computes the prediction polynomial, A, based on the
%   reflection coefficients, K.
%
%   [A,Efinal] = RC2POLY(K,R0) returns the final prediction error, Efinal,
%   based on the zero lag autocorrelation, R0.
%
%   % Example:
%   %   Consider a lattice IIR filter given by reflection coefficients,
%   %   k = [0.3090    0.9800    0.0031    0.0082   -0.0082], and give its
%   %   equivalent prediction filter representation.
%
%   k = [0.3090    0.9800    0.0031    0.0082   -0.0082];
%   a = rc2poly(k)      % Gives prediction polynomial
%
%   See also POLY2RC, RC2AC, AC2RC, AC2POLY, POLY2AC.

%   References: S. Kay, Modern Spectral Estimation,
%               Prentice Hall, N.J., 1987, Chapter 6.

%   Copyright 1988-2013 The MathWorks, Inc.

if (size(kr,1) > 1) && (size(kr,2) > 1)
    error(message('signal:rc2poly:inputnotsupported'));
end

% Initialize the recursion
kr = kr(:);               % Force kr to be a column vector.

p = length(kr);           % p is the order of the prediction polynomial.
a = [1 kr(1)];            % a is a true polynomial.

if (nargout == 2) && (nargin < 2),
    error(message('signal:rc2poly:SignalErr'));
end

% At this point nargin will be either 1 or 2
if nargin < 2,
    e0 = 0;  % Default value when e0 is not specified
else
    e0 = R0;
end

% Cast to enforce Precision rules
% This function accepts non-float numeric types but there are no enforced
% rules for the arithmetic in those cases. 
if any([signal.internal.sigcheckfloattype(kr,'single','rc2poly','K','allownumeric')...
    signal.internal.sigcheckfloattype(e0,'single','rc2poly','R0','allownumeric')])
  kr = single(kr);
  e0 = single(e0);
end

e(1) = e0.*(1-kr(1)'.*kr(1));

% Continue the recursion for k=2,3,...,p, where p is the order of the
% prediction polynomial.

for k = 2:p,
    [a,e(k)] = levup(a,kr(k),e(k-1)); %#ok<AGROW>
end

efinal = e(end);

% [EOF] rc2poly.m

