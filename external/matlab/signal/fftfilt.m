function y = fftfilt(b,x,nfft)
%FFTFILT Overlap-add method for FIR filtering using FFT.
%   Y = FFTFILT(B,X) filters X, with the FIR filter specified by the vector
%   of coefficients B, using the overlap/add method, and internal
%   parameters (FFT size and block length) that guarantee efficient
%   execution.
%   
%   Y = FFTFILT(B,X,N) allows you to have some control over the internal
%   parameters, by using an FFT of at least N points.
%
%   Y = FFTFILT(D,X,...) filters X with the FIR digital filter D. You
%   design a digital filter, D, by calling the <a href="matlab:help designfilt">designfilt</a> function.
%
%   If X is a matrix, FFTFILT filters its columns.  If B is a matrix,
%   FFTFILT applies the filter in each column of B to the signal vector X.
%   If B and X are both matrices with the same number of columns, then the
%   i-th column of B is used to filter the i-th column of X.
%
%   It is advantageous to use FFTFILT instead of FILTER when the signal is
%   relatively large.  FILTER performs N multiplications for each sample in
%   X where N is the filter length.  FFTFILT performs 2 FFT operations at
%   the cost of L*log2(L)/2 where L is the block length.  It then performs
%   L pointwise multiplications for a total cost of L*(1+log2(L))
%   multiplications.  The cost ratio is therefore L*(1+log2(L))/(N*L) =>
%   (1+log2(L))/N which is approximately log2(L)/N.  Therefore FFTFILT
%   becomes advantageous when log2(L) is less than N.
%
%   % Example 1:
%   %   Construct a Signal and filter it with a 10 point averaging filter
%   %   using fftfilt. 
%
%   fs = 100;                               % Sampling frequency
%   t = 0:1/fs:1;                           % Time vector
%   x = sin(2*pi*t*3)+.25*sin(2*pi*t*40);   % Input Signal
%   b = ones(1,10)/10;  % 10 point averaging filter
%   y = fftfilt(b,x);   % FIR filtering using overlap-add method
%   plot(t,x,t,y,'--');
%   legend('Original Signal','Filtered Signal')
%
%   % Example 2:
%   %   Use the designfilt function to design a lowpass FIR digital filter 
%   %   with order 350 and cutoff frequency of 150 Hz. The sample rate is 
%   %   1.5 KHz. Filter a long vector of data using the overlap-add method 
%   %   to increase speed.
%
%   D = designfilt('lowpassfir', 'FilterOrder', 350, ...
%    'CutoffFrequency', 150, 'SampleRate', 1500);
%
%   data = randn(10e6,1);  
%   y = fftfilt(D,data);
%
%   See also FILTER, FILTFILT.

%   --- Algorithmic details ---
%   The overlap/add algorithm convolves B with blocks of X, and adds
%   the overlapping output blocks.  It uses the FFT to compute the
%   convolution.
% 
%   Particularly for long FIR filters and long signals, this algorithm is 
%   MUCH faster than the equivalent numeric function FILTER(B,1,X).
%
%   Y = FFTFILT(B,X) -- If you leave N unspecified:   (RECOMMENDED)
%       Usually, length(X) > length(B).  Here, FFTFILT chooses an FFT 
%       length (N) and block length (L) which minimize the number of 
%       flops required for a length-N FFT times the number of blocks
%       ceil(length(X)/L).  
%       If length(X) <= length(B), FFTFILT uses a single FFT of length
%       nfft = 2^nextpow2(length(B)+length(X)-1), essentially computing 
%       ifft(fft(B,nfft).*fft(X,nfft)).
%
%   Y = FFTFILT(B,X,N) -- If you specify N:
%       In this case, N must be at least length(B); if it isn't, FFTFILT 
%       sets N to length(B).  Then, FFTFILT uses an FFT of length 
%       nfft = 2^nextpow2(N), and block length L = nfft - length(B) + 1. 
%       CAUTION: this can be VERY inefficient, if L ends up being small.

%   Author(s): L. Shure, 7-27-88
%              L. Shure, 4-25-90, revised
%              T. Krauss, 1-14-94, revised
%   Copyright 1988-2013 The MathWorks, Inc.

%   Reference:
%      A.V. Oppenheim and R.W. Schafer, Digital Signal 
%      Processing, Prentice-Hall, 1975.

narginchk(2,3);

m = size(x, 1);
if m == 1
    x = x(:);    % turn row into a column
end

nx = size(x,1);

if min(size(b))>1
   if (size(b,2)~=size(x,2))&&(size(x,2)>1)
      error(message('signal:fftfilt:InvalidDimensions'))
   end
else
   b = b(:);   % make input a column
end
nb = size(b,1);

if nargin < 3
% figure out which nfft and L to use
    if nb >= nx || nb > 2^20    % take a single FFT in this case
        nfft = 2^nextpow2(nb+nx-1);
        L = nx;
    else
        fftflops = [ 18 59 138 303 660 1441 3150 6875 14952 32373 69762 ...
       149647 319644 680105 1441974 3047619 6422736 13500637 28311786 59244791];
        n = 2.^(1:20); 
        validset = find(n>(nb-1));   % must have nfft > (nb-1)
        n = n(validset); 
        fftflops = fftflops(validset);
        % minimize (number of blocks) * (number of flops per fft)
        L = n - (nb - 1);
        [dum,ind] = min( ceil(nx./L) .* fftflops ); %#ok
        nfft = n(ind);
        L = L(ind);
    end

else  % nfft is given
  % Cast to enforce precision rules
    nfft = signal.internal.sigcasttofloat(nfft,'double','fftfilt','N','allownumeric');
    if nfft < nb
        nfft = nb;
    end
    nfft = 2.^(ceil(log(nfft)/log(2))); % force this to a power of 2 for speed
    L = nfft - nb + 1;
end

% Check the input data type. Single precision is not supported.
try
    chkinputdatatype(b,x,nfft);
catch ME
    throwAsCaller(ME);
end

B = fft(b,nfft);
if length(b)==1,
     B = B(:);  % make sure fft of B is a column (might be a row if b is scalar)
end
if size(b,2)==1
    B = B(:,ones(1,size(x,2)));  % replicate the column B 
end
if size(x,2)==1
    x = x(:,ones(1,size(b,2)));  % replicate the column x 
end
y = zeros(size(x));

istart = 1;
while istart <= nx
    iend = min(istart+L-1,nx);
    if (iend - istart) == 0
        X = x(istart(ones(nfft,1)),:);  % need to fft a scalar
    else
        X = fft(x(istart:iend,:),nfft);
    end
    Y = ifft(X.*B);
    yend = min(nx,istart+nfft-1);
    y(istart:yend,:) = y(istart:yend,:) + Y(1:(yend-istart+1),:);
    istart = istart + L;
end

if ~(any(imag(b(:))) || any(imag(x(:))))
	y = real(y);
end

if (m == 1)&&(size(y,2) == 1)
    y = y(:).';    % turn column back into a row
end

