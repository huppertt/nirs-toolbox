function [a,efinal]=ac2poly(R)
%AC2POLY  Convert autocorrelation sequence to prediction polynomial.
%   [A,Efinal]=AC2POLY(R) returns the prediction polynomial, A, and the final 
%   prediction error, Efinal, based on the autocorrelation sequence, R.
%
%   % Example:
%   %   Consider a autocorrelation sequence and find the corresponding 
%   %   prediction filter polynomial.
%
%   r = [5.0000 -1.5450 -3.9547 3.9331 1.4681 -4.7500];
%   [a,efinal] = ac2poly(r)             % Prediction polynomial
% 
%   See also POLY2AC, POLY2RC, RC2POLY, RC2AC, AC2RC. 


%   References: S. Kay, Modern Spectral Estimation,
%               Prentice Hall, N.J., 1987, Chapter 6.
%
%   Author(s): A. Ramasubramanian
%   Copyright 1988-2002 The MathWorks, Inc.

% Check the input data type. Single precision is not supported.
try
    chkinputdatatype(R);
catch ME
    throwAsCaller(ME);
end

% Use levinson recursion for this. Note that matrix inversion 
% could be used for small orders.

[a,efinal] = levinson(R);

% [EOF] ac2poly.m         


