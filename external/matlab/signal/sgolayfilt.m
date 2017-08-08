function y=sgolayfilt(x,k,F,W,DIM)
%SGOLAYFILT Savitzky-Golay Filtering.
%   SGOLAYFILT(X,K,F) smoothes the signal X using a Savitzky-Golay 
%   (polynomial) smoothing filter.  The polynomial order, K, must
%   be less than the frame size, F, and F must be odd.  The length 
%   of the input X must be >= F.  If X is a matrix, the filtering
%   is done on the columns of X.
%
%   Note that if the polynomial order K equals F-1, no smoothing
%   will occur.
%
%   SGOLAYFILT(X,K,F,W) specifies a weighting vector W with length F
%   containing real, positive valued weights employed during the
%   least-squares minimization. If not specified, or if specified as
%   empty, W defaults to an identity matrix.
%
%   SGOLAYFILT(X,K,F,[],DIM) or SGOLAYFILT(X,K,F,W,DIM) operates along
%   the dimension DIM.
%
%   % Example:
%   %   Smooth the mtlb signal by applying a cubic Savitzky-Golay filter 
%   %   to data frames of length 41.
%
%   load mtlb                        % Load data
%   smtlb = sgolayfilt(mtlb,3,41);   % Apply 3rd-order filter
%   h1 = subplot(2,1,1)
%   plot([1:2000],mtlb(1:2000)); axis([0 2000 -4 4]);
%   title('Original Data'); grid;
%   h2 = subplot(2,1,2)
%   plot([1:2000],smtlb(1:2000)); axis([0 2000 -4 4]);
%   title('Filtered Data'); grid;
%   linkaxes([h1,h2],'x')         % Linking subplots
%
%   See also SGOLAY, MEDFILT1, FILTER

%   References:
%     [1] Sophocles J. Orfanidis, INTRODUCTION TO SIGNAL PROCESSING,
%              Prentice-Hall, 1995, Chapter 8.

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(3,5,nargin,'struct'));

% Check if the input arguments are valid
if round(F) ~= F, error(message('signal:sgolayfilt:MustBeIntegerFrameLength')), end
if rem(F,2) ~= 1, error(message('signal:sgolayfilt:SignalErr')), end
if round(k) ~= k, error(message('signal:sgolayfilt:MustBeIntegerPolyDegree')), end
if k > F-1, error(message('signal:sgolayfilt:InvalidRangeDegree')), end

if nargin < 4 | isempty(W), %#ok
   % No weighting matrix, make W an identity
   W = ones(F,1);
else
   % Check for right length of W
   if length(W) ~= F, error(message('signal:sgolayfilt:InvalidDimensionsWeight')),end
   % Check to see if all elements are positive
   if min(W) <= 0, error(message('signal:sgolayfilt:InvalidRangeWeight')), end
end

if nargin < 5, DIM = []; end

% Check the input data type. Single precision is not supported.
try
    chkinputdatatype(x,k,F,W,DIM);
catch ME
    throwAsCaller(ME);
end

% Compute the projection matrix B
B = sgolay(k,F,W);

if ~isempty(DIM) && DIM > ndims(x)
	error(message('signal:sgolayfilt:InvalidDimensionsInput', 'X'))
end

% Reshape X into the right dimension.
if isempty(DIM)
	% Work along the first non-singleton dimension
	[x, nshifts] = shiftdim(x);
else
	% Put DIM in the first dimension (this matches the order 
	% that the built-in filter function uses)
	perm = [DIM,1:DIM-1,DIM+1:ndims(x)];
	x = permute(x,perm);
end

if size(x,1) < F, error(message('signal:sgolayfilt:InvalidDimensionsTooSmall')), end

% Preallocate output
y = zeros(size(x));

% Compute the transient on
y(1:(F+1)/2-1,:) = flipud(B((F-1)/2+2:end,:))*flipud(x(1:F,:));

% Compute the steady state output
ytemp = filter(B((F-1)./2+1,:),1,x);
y((F+1)/2:end-(F+1)/2+1,:) = ytemp(F:end,:);

% Compute the transient off
y(end-(F+1)/2+2:end,:) = flipud(B(1:(F-1)/2,:))*flipud(x(end-(F-1):end,:));

% Convert Y to the original shape of X
if isempty(DIM)
	y = shiftdim(y, -nshifts);
else
	y = ipermute(y,perm);
end

