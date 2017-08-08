function Y = goertzel(X,INDVEC,DIM)
%GOERTZEL Second-order Goertzel algorithm.  
%   GOERTZEL(X,INDVEC) computes the discrete Fourier transform (DFT)
%   of X at indices contained in the vector INDVEC, using the
%   second-order Goertzel algorithm.  The indices must be integer values
%   from 1 to N where N is the length of the first non-singleton dimension.
%   If empty or omitted, INDVEC is assumed to be 1:N.
%
%   GOERTZEL(X,[],DIM) or GOERTZEL(X,INDVEC,DIM) computes the DFT along 
%   the dimension DIM.
%
%   In general, GOERTZEL is slower than FFT when computing all the possible
%   DFT indices, but is most useful when X is a long vector and the DFT 
%   computation is required for only a subset of indices less than
%   log2(length(X)).  Indices 1:length(X) correspond to the frequency span
%   [0, 2*pi) radians.
%
%   EXAMPLE:
%      % Resolve the 1.24 kHz and 1.26 kHz components in the following
%      % noisy cosine which also has a 10 kHz component.
%      Fs = 32e3;   t = 0:1/Fs:2.96;
%      x  = cos(2*pi*t*10e3)+cos(2*pi*t*1.24e3)+cos(2*pi*t*1.26e3)...
%           + randn(size(t));
%
%      N = (length(x)+1)/2;
%      f = (Fs/2)/N*(0:N-1);              % Generate frequency vector
%      indxs = find(f>1.2e3 & f<1.3e3);   % Find frequencies of interest
%      X = goertzel(x,indxs);
%      
%      plot(f(indxs)/1e3,20*log10(abs(X)/length(X)));
%      title('Mean Squared Spectrum');
%      xlabel('Frequency (kHz)');
%      ylabel('Power (dB)');
%      grid on;
%      set(gca,'XLim',[f(indxs(1)) f(indxs(end))]/1e3);
%
%   See also FFT, FFT2.

%   Copyright 1988-2013 The MathWorks, Inc.

%   Reference:
%     C.S.Burrus and T.W.Parks, DFT/FFT and Convolution Algorithms, 
%     John Wiley & Sons, 1985

if ~(exist('goertzelmex', 'file') == 3), 
	error(message('signal:goertzel:NotSupported'));
end

narginchk(1,3);
if nargin < 2, INDVEC = []; end
if nargin < 3, DIM = []; end

% Cache input data type to cast output at end of computations
isDataSingle = false;
if isa(X,'single')
  isDataSingle = true;
end
X = signal.internal.sigcasttofloat(X,'double','goertzel','X');
INDVEC = signal.internal.sigcasttofloat(INDVEC,'double','goertzel',...
  'INDVEC','allownumeric');
DIM = signal.internal.sigcasttofloat(DIM,'double','goertzel','DIM',...
  'allownumeric');

if ~isempty(DIM) && DIM > ndims(X)
	error(message('signal:goertzel:InvalidDimensions'))
end

% Reshape X into the right dimension.
if isempty(DIM)
	% Work along the first non-singleton dimension
	[X, nshifts] = shiftdim(X);
else
	% Put DIM in the first dimension (this matches the order 
	% that the built-in filter function uses)
	perm = [DIM,1:DIM-1,DIM+1:ndims(X)];
	X = permute(X,perm);
end

% Verify that the indices in INDVEC are valid.
siz = size(X);
if isempty(INDVEC),
	INDVEC = 1:siz(1); % siz(1) is the number of rows of X
else
	INDVEC = INDVEC(:);
	if max(INDVEC) > siz(1),
		error(message('signal:goertzel:IdxGtBound'));
	elseif min(INDVEC) < 1
		error(message('signal:goertzel:IdxLtBound'));
	elseif all(INDVEC-fix(INDVEC))
		error(message('signal:goertzel:MustBeInteger'));
	end
end

% Initialize Y with the correct dimension
Y = zeros([length(INDVEC),siz(2:end)]); 

% Call goertzelmex 
for k = 1:prod(siz(2:end)),
	Y(:,k) = goertzelmex(X(:,k),INDVEC-1);
end

% Convert Y to the original shape of X
if isempty(DIM)
	Y = shiftdim(Y, -nshifts);
else
	Y = ipermute(Y,perm);
end

% If X is single, then relevant output should be single
if isDataSingle
  Y = single(Y);
end

% [EOF] goertzel.m
