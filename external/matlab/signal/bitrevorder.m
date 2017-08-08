function [y,idx] = bitrevorder(x)
%BITREVORDER Permute input into bit-reversed order.
%   Y = BITREVORDER(X) returns an array Y of the same size as X, with its
%   first non-singleton dimension permuted to be in bit-reversed order.
%
%   [Y,I] = BITREVORDER(X) returns the bit-reversed index vector I, such
%   that, for a vector Y, Y=X(I).
%
%   This operation is useful to pre-order a vector of filter coefficients
%   for use in frequency-domain filtering algorithms, in which the FFT and
%   IFFT transforms are computed without bit-reversed ordering for improved
%   run-time efficiency.
%
%   Calling BITREVORDER(X) is equivalent to calling DIGITREVORDER(X, 2).
%
%   EXAMPLE:
%       x = (0:15).';
%       y = bitrevorder(x); % radix-2 (bit-reversed) ordering of x
%       [x y]
%
%   See also DIGITREVORDER, FFT, IFFT.

% Copyright 1988-2010 The MathWorks, Inc.

[y, idx] = digitrevorder(x, 2);
