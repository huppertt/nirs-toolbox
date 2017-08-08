function c = cconv(a,b,N)
%CCONV Modulo-N circular convolution.
%   C = CCONV(A, B, N) circularly convolves vectors A and B.  The resulting
%   vector is length N. If omitted, N defaults to LENGTH(A)+LENGTH(B)-1.
%   When N = LENGTH(A)+LENGTH(B)-1, the circular convolution is equivalent
%   to the linear convolution computed by CONV.
%   
%   % Example #1: Mod-4 circular convolution
%   a = [2 1 2 1];
%   b = [1 2 3 4];
%   c = cconv(a,b,4)  
%
%   % Example #2: Circular convolution as a fast linear convolution
%   a = [1 2 -1 1];
%   b = [1 1 2 1 2 2 1 1];
%   c = cconv(a,b,11)  
%   cref = conv(a,b)
%
%   % Example #3: Circular cross-correlation
%   a = [1 2 2 1]+1i;
%   b = [1 3 4 1]-2*1i;
%   c = cconv(a,conj(fliplr(b)),7)
%   cref = xcorr(a,b)
%   
%   See also CONV, XCORR

%   Copyright 2006-2010 The MathWorks, Inc.

%   Reference:
%     Sophocles J. Orfanidis, Introduction to Signal Processing, 
%     Prentice-Hall, 1996

error(nargchk(2,3,nargin,'struct'))

na = length(a);
nb = length(b);

if na ~= numel(a) || nb ~= numel(b)
  error(message('signal:cconv:AorBNotVector'));
end


if nargin<3,
    N = na+nb-1;
end
c = ifft(fft(datawrap(a,N),N).*fft(datawrap(b,N),N));

% [EOF]
