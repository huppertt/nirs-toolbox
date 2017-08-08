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
%   See also CCONV, GPUARRAY

%   Copyright 2012 The MathWorks, Inc.

%   Reference:
%     Sophocles J. Orfanidis, Introduction to Signal Processing, 
%     Prentice-Hall, 1996

narginchk(2,3);

if ~ (isa(a, 'gpuArray') || isa(b, 'gpuArray') ),
  %if only N is a gpuArray, dispatch to cpu version
  N = gather(N);
  c = cconv(a,b,N);
  return;
else
  %put data in the right place
  if nargin > 2,
    N = gather(N);
  end
  a = gpuArray(a);
  b = gpuArray(b);
end

na = length(a);
nb = length(b);

if na ~= numel(a) || nb ~= numel(b)
  error(message('signal:cconv:AorBNotVector'));
end

if isreal(a) && isreal(b)
    symFlag = 'symmetric';
else
    symFlag = 'nonsymmetric';
end

if nargin<3,
    N = na+nb-1;
end
c = ifft(fft(privDataWrap(a,N),N).*fft(privDataWrap(b,N),N), symFlag);


function x = privDataWrap(x,nfft)
nx = size(x, 2);
% Reshape into multiple columns (data segments) of length nfft.
% If insufficient data points are available, zeros are appended.
% Sum across the columns (data segments).
x = sum(privBuffer(x,nfft),2);
% Reshape vector as necessary:
if (nx~=1), 
    x=x.'; 
end

function y = privBuffer(x,n)

rem = mod(numel(x),n);

if (rem ~= 0),
    pad = n-rem;
    y = [reshape(x,[],1); ...
        gpuArray.zeros(pad,1)];
    y = reshape(y,n, []);
else
    y = reshape(x,n, []);
end

% [EOF]
