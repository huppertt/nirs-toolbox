function D = dftmtx(n)
%DFTMTX Discrete Fourier transform matrix.
%   DFTMTX(N) is the N-by-N complex matrix of values around
%   the unit-circle whose inner product with a column vector
%   of length N yields the discrete Fourier transform of the
%   vector.  If X is a column vector of length N, then
%   DFTMTX(N)*X yields the same result as FFT(X); however, 
%   FFT(X) is more efficient.
%
%   The inverse discrete Fourier transform matrix is
%   CONJ(DFTMTX(N))/N.
%
%   An interesting example is 
%
%     D = dftmtx(4)
%
%   which returns
%
%     D = [1   1   1   1
%          1  -i  -1   i     
%          1  -1   1  -1
%          1   i  -1  -i]
%
%   which illustrates why no multiplications are necessary for
%   4-point DFTs.
%
%   See also FFT and IFFT.

%   References:
%     [1] Alan V. Oppenheim and Ronald W. Schafer, Discrete-Time
%     Signal Processing, Prentice Hall, 1989.
%
%     [2] Charles Van Loan, Computational Frameworks for the Fast
%     Fourier Transform, SIAM, 1992.

%   Copyright 1988-2013 The MathWorks, Inc.

n = signal.internal.sigcasttofloat(n,'double','dftmtx','N',...
  'allownumeric');

D = fft(eye(n));

% Explicitly, the calculation is:
%  f = 2*pi/n;                 % Angular increment.
%  w = (0:f:2*pi-f/2).' * 1i;  % Column.
%  x = 0:n-1;                  % Row.
%  D = exp(-w*x);              % Exponentiation of outer product.
