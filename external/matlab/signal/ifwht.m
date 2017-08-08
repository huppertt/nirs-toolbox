function y = ifwht(varargin)
%IFWHT Fast Inverse Discrete Walsh-Hadamard Transform
%   Y = IFWHT(X) returns the inverse discrete Walsh-Hadamard transform of
%   X. The inverse transform values are stored in Y. If X is a matrix, the 
%   function operates on each column.
%
%   Y = IFWHT(X,N) is the N-point inverse Walsh-Hadamard transform of X
%   where N must be a power of two. X is padded with zeros if it has less 
%   than N points and truncated if it has more. The default value of N is 
%   equal to the length of the vector X if it is a power of two or the next
%   power of two greater than the signal vector length. The function
%   returns an error if N is not a power of 2.
%
%   Y = IDWHT(X,N,ORDERING) is the N-point inverse discrete Walsh-Hadamard
%   transform of the vector X where the output samples values are output in
%   the desired order as specified by ORDERING. ORDERING can be 'sequency',
%   'hadamard' or 'dyadic'. The default ORDERING is 'sequency'.
%    
%   EXAMPLE:
%           x = rand(16);
%           y = fwht(x);
%           xHat = ifwht(y); % Inverse transformation should reproduce x
%
%   For more information see the <a href="matlab:web([matlabroot,'\toolbox\signal\sigdemos\html\walshhadamarddemo.html'])">Walsh-Hadamard Transform Demo</a> or enter "doc ifwht"
%   at the MATLAB command line.
%
%   See also FWHT, FFT, IFFT, DCT, IDCT, DWT, IDWT.

%   Copyright 2008-2013 The MathWorks, Inc.

narginchk(1,3)
% Since the forward and inverse transforms are exactly identical
% operations, FWHT is used to perform inverse transform
y = fwht(varargin{:});
% Perform scaling 
[m,n] = size(y);
if m == 1 % column vector
    m = n;
end
y = y .* m;
