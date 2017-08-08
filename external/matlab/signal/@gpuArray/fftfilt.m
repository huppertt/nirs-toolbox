function y = fftfilt(b,x,nfft)
%FFTFILT Overlap-add method of FIR filtering using FFT.
%   Y = FFTFILT(B,X) filters X with the FIR filter B using the overlap/add
%   method, using internal parameters (FFT size and block length) which
%   guarantee efficient execution.
%   
%   Y = FFTFILT(B,X,N) allows you to have some control over the internal
%   parameters, by using an FFT of at least N points.
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
%   See also FILTFILT, GPUARRAY.

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

%   Copyright 2012 The MathWorks, Inc.

%   Reference:
%      A.V. Oppenheim and R.W. Schafer, Digital Signal 
%      Processing, Prentice-Hall, 1975.

narginchk(2,3); %check input argument count

if ~ (isa(x, 'gpuArray') || isa(b, 'gpuArray') ),
  %if only nfft is a gpuArray, dispatch to cpu version
  nfft = gather(nfft);
  y = fftfilt(b,x,nfft);
  return;
else
  %put data in the right place
  if nargin > 2,
    nfft = gather(nfft);
  end
  x = gpuArray(x);
  b = gpuArray(b);
end

m = size(x, 1);
if m == 1
    x = reshape(x,[],1);    % turn row into a column
end

nx = size(x,1);

if min(size(b))>1
   if (size(b,2)~=size(x,2))&&(size(x,2)>1)
      error(message('signal:fftfilt:InvalidDimensions'));
   end
else
   b = reshape(b, [],1);   % make input a column
end
nb = size(b,1);

% This GPU algorithm handles empty arrays differently than the CPU.
% To avoid differences, explicitly handle the empties here and return.
if ( isempty(x) || isempty(b) )
    y = gpuArray.zeros(size(x));
    return;
end



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
    if nfft < nb
        nfft = nb;
    end
    nfft = 2.^(ceil(log(nfft)/log(2))); % force this to a power of 2 for speed
    L = nfft - nb + 1;
end



% Check the input data type. Single precision is not supported.       
if ~(checkDoubles(b) && checkDoubles(x) && checkDoubles(nfft))
    error(message('signal:gpu:OnlyDoubles'));
end

if isreal(b) && isreal(x)
    ifftflag  = 'symmetric';
else
    ifftflag = 'nonsymmetric';
end


B = fft(b,nfft,1);


%Buffer x into sections of length L by inserting zeros and pad with an
%extra column of zeros. FFT this buffered x, multiply in frequency
%domain and then IFFT. The bottom L:nfft rows are the overlap. Using
%indexing, add the overlap from the bottom of one column to the top of the
%next. Then remove the overlap rows. Reshape to get back to the proper
%output size.


xbuf = privBuffer(x, L);
X = fft(xbuf, nfft,1);
if length(x)==1,
    X = X(:);
end

[Xr, ~] = size(X);
colsOut = max(size(b,2), size(x,2));

X = reshape(X, Xr, [], size(x,2));
B = reshape(B, size(B,1), 1, size(b,2));
Y = bsxfun(@times, X,B);
Y = reshape(Y, size(Y,1), []);

y = ifft(Y, nfft, ifftflag);

%overlap = y((L+1):nfft, 1:(size(y,2)-1));
overlap_subsr = substruct( '()', { ((L+1):nfft ), (1:(size(y,2)-1))} );
overlap = subsref(y, overlap_subsr);

%tmp = y(1:size(overlap,1), 2:size(y,2)) + overlap;
y_subs = substruct( '()', { 1:size(overlap,1), ...
                            2:size(y,2) } );
tmp = subsref(y, y_subs) + overlap;


%y(1:size(overlap,1), 2:size(y,2)) = tmp;
y = subsasgn(y, y_subs, tmp);

%now remove overlap
%y = y(1:L, 1:size(y,2))
ygood_subs = substruct( '()', {1:L, 1:size(y,2)});
y = subsref(y, ygood_subs);



y = reshape(y, [], colsOut);

%remove tail
%y((size(x,1)+1):size(y,1), :) = []
ynotail_subs = substruct( '()', {(size(x,1)+1):size(y,1), ':' } );
y = subsasgn(y, ynotail_subs, []);



if (m == 1)&&(size(y,2) == 1)
    y = reshape(y,1, []);    % turn column back into a row
end

end





function y = privBuffer(x,n)
% Add a matrix of zeros to the bottom of x,
% then n more zeros per column. Reshape to the appropriate size.
% At the end the matrix y should look like:
%   y = [ x(a:b,1)  x(b+1:c, 1) ... 0  x(a:b,2) x(b+1:c,2) ...0
%            0           0       ...0      0        0      ...0
%                   ....
%             0          0          0      0        0         0
%

  rem = mod(size(x,1),n);

  if (rem ~= 0),
      pad = n-rem;
      y = [x ;...
          gpuArray.zeros(pad,size(x,2))];
  else
      y = x;
  end

  y = [y; gpuArray.zeros(n, size(x,2))];  %one more column worth so the 
                                   %overlap doesn't cross from one frame to
                                   %the next.
  y = reshape(y, n, []);

end

function tf = checkDoubles(z)
    if isa(z, 'double')
        tf = true;
    elseif isa(z, 'gpuArray')
        tf = strcmpi(classUnderlying(z), 'double');
    end
end
        
