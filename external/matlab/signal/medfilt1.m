function y = medfilt1(x,n,blksz,DIM)
%MEDFILT1  One dimensional median filter.
%   Y = MEDFILT1(X,N) returns the output of the order N, one dimensional
%   median filtering of X.  Y is the same size as X; for the edge points,
%   zeros are assumed to the left and right of X.  If X is a matrix,
%   then MEDFILT1 operates along the columns of X.
%
%   If you do not specify N, MEDFILT1 uses a default of N = 3.
%   For N odd, Y(k) is the median of X( k-(N-1)/2 : k+(N-1)/2 ).
%   For N even, Y(k) is the median of X( k-N/2 : k+N/2-1 ).
%
%   Y = MEDFILT1(X,N,BLKSZ) uses a for-loop to compute BLKSZ ("block size") 
%   output samples at a time.  Use this option with BLKSZ << LENGTH(X) if 
%   you are low on memory (MEDFILT1 uses a working matrix of size
%   N x BLKSZ).  By default, BLKSZ == LENGTH(X); this is the fastest
%   execution if you have the memory for it.
%
%   For matrices and N-D arrays, Y = MEDFILT1(X,N,[],DIM) or 
%   Y = MEDFILT1(X,N,BLKSZ,DIM) operates along the dimension DIM.
%   
%   % Example:
%   %   Construct a noisy signal and apply a 10th order one-dimensional 
%   %   median filter to it.
%
%   fs = 100;                               % Sampling rate                                   
%   t = 0:1/fs:1;                           % Time vector
%   x = sin(2*pi*t*3)+.25*sin(2*pi*t*40);   % Noise Signal - Input
%   y = medfilt1(x,10);                     % Median filtering - Output
%   plot(t,x,'k',t,y,'r'); grid;            % Plot 
%   legend('Original Signal','Filtered Signal')
%
%   See also MEDFILT2, MEDIAN, FILTER, SGOLAYFILT.
%
%   Note:  MEDFILT2 is in the Image Processing Toolbox.

%   Author(s): L. Shure and T. Krauss, 8-3-93
%   Copyright 1988-2014 The MathWorks, Inc.

% Validate number of input arguments
error(nargchk(1,4,nargin,'struct'));
if nargin < 2, n = []; end
if nargin < 3, blksz = []; end
if nargin < 4, DIM = []; end

% Check the input data type. Single precision is not supported.
try
    chkinputdatatype(x,n,blksz,DIM);
catch ME
    throwAsCaller(ME);
end

% Check if the input arguments are valid
if isempty(n)
  n = 3;
end

if ~isempty(DIM) && DIM > ndims(x)
	error(message('signal:medfilt1:InvalidDimensions'))
end

% Reshape x into the right dimension.
if isempty(DIM)
	% Work along the first non-singleton dimension
	[x, nshifts] = shiftdim(x);
else
	% Put DIM in the first (row) dimension (this matches the order 
	% that the built-in filter function uses)
	perm = [DIM,1:DIM-1,DIM+1:ndims(x)];
	x = permute(x,perm);
end

% Verify that the block size is valid.
siz = size(x);
if isempty(blksz),
	blksz = siz(1); % siz(1) is the number of rows of x (default)
else
	blksz = blksz(:);
end

% Initialize y with the correct dimension
y = zeros(siz); 

% Call medfilt1D (vector)
for i = 1:prod(siz(2:end)),
	y(:,i) = medfilt1D(x(:,i),n,blksz);
end

% Convert y to the original shape of x
if isempty(DIM)
	y = shiftdim(y, -nshifts);
else
	y = ipermute(y,perm);
end


%-------------------------------------------------------------------
%                       Local Function
%-------------------------------------------------------------------
function y = medfilt1D(x,n,blksz)
%MEDFILT1D  One dimensional median filter.
%
% Inputs:
%   x     - vector
%   n     - order of the filter
%   blksz - block size

nx = length(x);
if rem(n,2)~=1    % n even
    m = n/2;
else
    m = (n-1)/2;
end
X = [zeros(m,1); x; zeros(m,1)];
y = zeros(nx,1);

% Work in chunks to save memory
indr = (0:n-1)';
indc = 1:nx;
for i=1:blksz:nx
    ind = indc(ones(1,n),i:min(i+blksz-1,nx)) + ...
          indr(:,ones(1,min(i+blksz-1,nx)-i+1));
    xx = reshape(X(ind),n,min(i+blksz-1,nx)-i+1);
    y(i:min(i+blksz-1,nx)) = median(xx,1);
end
