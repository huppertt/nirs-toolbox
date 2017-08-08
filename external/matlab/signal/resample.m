function  varargout = resample(varargin)
%RESAMPLE  Resample uniform or nonuniform data to a new fixed rate.
%   Y = RESAMPLE(X,P,Q) resamples the values, X, of a uniformly sampled
%   signal at P/Q times the original sample rate using a polyphase
%   anti-aliasing filter. If X is a matrix, then RESAMPLE treats each
%   column as an independent channel.
%
%   In its filtering process, RESAMPLE assumes the samples at times before
%   and after the given samples in X are equal to zero. Thus large
%   deviations from zero at the end points of the sequence X can cause
%   inaccuracies in Y at its end points.
%
%   [Y,Ty] = RESAMPLE(X,Tx) resamples the values, X, of a signal sampled at
%   the instants specified in vector Tx. RESAMPLE interpolates X linearly
%   onto a vector of uniformly spaced instants, Ty, with the same endpoints
%   and number of samples as Tx.  NaNs are treated as missing data and are
%   ignored.
%
%   [Y,Ty] = RESAMPLE(X,Tx,Fs) uses interpolation and an anti-aliasing
%   filter to resample the signal at a uniform sample rate, Fs.
% 
%   [Y,Ty] = RESAMPLE(X,Tx,Fs,P,Q) interpolates X to an intermediate
%   uniform grid with sample rate equal Q*Fs/P and filters the result
%   using UPFIRDN to upsample by P and downsample by Q.  Specify P and Q
%   so that Q*Fs/P is least twice as large as the highest frequency in the
%   input signal.  
%
%   [Y,Ty] = RESAMPLE(X,Tx,...,METHOD) specifies the interpolation method.
%   The default is linear interpolation.  Available methods are:
%     'linear' - linear interpolation
%     'pchip'  - shape-preserving piecewise cubic interpolation
%     'spline' - piecewise cubic spline interpolation
%
%   Y = RESAMPLE(...,P,Q,N) uses a weighted sum of 2*N*max(1,Q/P) samples
%   of X to compute each sample of Y.  The length of the FIR filter
%   RESAMPLE applies is proportional to N; by increasing N you will get
%   better accuracy at the expense of a longer computation time.  If you
%   don't specify N, RESAMPLE uses N = 10 by default.  If you let N = 0,
%   RESAMPLE performs a nearest neighbor interpolation; that is, the output
%   Y(n) is X(round((n-1)*Q/P)+1) ( Y(n) = 0 if round((n-1)*Q/P)+1 >
%   length(X) ).
%
%   Y = RESAMPLE(...,P,Q,N,BTA) uses BTA as the BETA design parameter for
%   the Kaiser window used to design the filter.  RESAMPLE uses BTA = 5 if
%   you don't specify a value.
%
%   Y = RESAMPLE(...,P,Q,B) uses B to filter X (after upsampling) if B is a
%   vector of filter coefficients.  RESAMPLE assumes B has odd length and
%   linear phase when compensating for the filter's delay; for even length
%   filters, the delay is overcompensated by 1/2 sample.  For non-linear
%   phase filters consider using UPFIRDN.
%
%   [Y,B] = RESAMPLE(X,P,Q,...) returns in B the coefficients of the filter
%   applied to X during the resampling process (after upsampling).
%
%   [Y,Ty,B] = RESAMPLE(X,Tx,...) returns in B the coefficients of the
%   filter applied to X during the resampling process (after upsampling).
%
%   % Example 1:
%   %   Resample a sinusoid at 3/2 the original rate.
%
%   tx = 0:3:300-3;         % Time vector for original signal
%   x = sin(2*pi*tx/300);   % Define a sinusoid 
%   ty = 0:2:300-2;         % Time vector for resampled signal        
%   y = resample(x,3,2);    % Change sampling rate
%   plot(tx,x,'+-',ty,y,'o:')
%   legend('original','resampled');
%   xlabel('Time')
%
%   % Example 2:
%   %   Resample a non-uniformly sampled sinusoid to a uniform 50 Hz rate.
%   
%   Fs = 50;
%   tx = linspace(0,1,21) + .012*rand(1,21);
%   x = sin(2*pi*tx);
%   [y, ty] = resample(x, tx, Fs);
%   plot(tx,x,'+-',ty,y,'o:')
%   legend('original','resampled');
%   xlabel('Time')
%   axis tight
%
%   See also UPFIRDN, INTERP, INTERP1, DECIMATE, FIRLS, KAISER, INTFILT.

%   NOTE-1: digital anti-alias filter is designed via windowing

%   Author(s): James McClellan, 6-11-93
%              Modified to use upfirdn, T. Krauss, 2-27-96
%   Copyright 1988-2014 The MathWorks, Inc.
%     

narginchk(2,8);

% fetch the users desired interpolation method (if any)
[method, varargin] = getInterpMethod(varargin);

if isscalar(varargin{2})
  % [...] = RESAMPLE(X,P,Q,...)
  % interpolation is not performed when using uniformly sampled data
  if ~strcmp(method,'')
    error(message('signal:resample:UnexpectedInterpolation',method));
  end
  [varargout{1:max(1,nargout)}] = uniformResample(varargin{:});
else
  % [...] = RESAMPLE(X,Tx,...)
  if strcmp(method,'')
    % use linear method by default
    method = 'linear';
  end
  [varargout{1:max(1,nargout)}] = nonUniformResample(method,varargin{:});
end


function [y, ty, h] = nonUniformResample(method,varargin)

% fetch non-NaN input samples in time-sorted order
[x, tx] = getSamples(varargin{1:2});
  
% obtain sample rate
if numel(varargin)>2
  fs = varargin{3};
  validateFs(fs);
else
  % compute the average sample rate if unspecified
  fs = (numel(tx)-1) / (tx(end)-tx(1));
end

% get rational sample ratio  
if numel(varargin)>4
  p = varargin{4};
  q = varargin{5};
  validateResampleRatio(p,q);
elseif numel(varargin)>3
  error(message('signal:resample:MissingQ'));
else
  % get a close rational approximation of the ratio between the desired
  % sampling rate and the average sample rate
  [p, q] = getResampleRatio(tx,fs);
end

% use cubic spline interpolation onto a uniformly sampled grid
% with a target sample rate of Q*Fs/P
tsGrid = p/(q*fs);
tGrid = tx(1):tsGrid:tx(end);


if isreal(x)
  xGrid = matInterp1(tx,x,tGrid,method);
else
  % compute real and imaginary channels independently
  realGrid = matInterp1(tx,real(x),tGrid,method);
  imagGrid = matInterp1(tx,imag(x),tGrid,method);
  xGrid = complex(realGrid,imagGrid);
end

% recover the desired sampling rate by resampling the grid at the 
% specified ratio
[y, h] = uniformResample(xGrid,p,q,varargin{6:end});

% create output time vector
if isvector(y)
  ty = tx(1) + (0:numel(y)-1)/fs;
else
  ty = tx(1) + (0:size(y,1)-1)/fs;
end  

% match dimensionality of output time vector to input
if iscolumn(tx)
  ty = ty(:);
end 

%-------------------------------------------------------------------------
function  [y, h] = uniformResample( x, p, q, N, bta )

if nargin < 5,  bta = 5;  end   %--- design parameter for Kaiser window LPF
if nargin < 4,   N = 10;   end

validateResampleRatio(p, q);

[p,q] = rat( p/q, 1e-12 );  %--- reduce to lowest terms 
   % (usually exact, sometimes not; loses at most 1 second every 10^12 seconds)
if (p==1) && (q==1)
    y = x; 
    h = 1;
    return
end
pqmax = max(p,q);
if length(N)>1      % use input filter
   L = length(N);
   h = N;
else                % design filter
   if( N>0 )
      fc = 1/2/pqmax;
      L = 2*N*pqmax + 1;
      h = firls( L-1, [0 2*fc 2*fc 1], [1 1 0 0]).*kaiser(L,bta)' ;
      h = p*h/sum(h);
   else
      L = p;
      h = ones(1,p);
   end
end

Lhalf = (L-1)/2;
isvect = any(size(x)==1);
if isvect
    Lx = length(x);
else
    Lx = size(x, 1);
end

% Need to delay output so that downsampling by q hits center tap of filter.
nz = floor(q-mod(Lhalf,q));
z = zeros(1,nz);
h = [z h(:).'];  % ensure that h is a row vector.
Lhalf = Lhalf + nz;

% Number of samples removed from beginning of output sequence 
% to compensate for delay of linear phase filter:
delay = floor(ceil(Lhalf)/q);

% Need to zero-pad so output length is exactly ceil(Lx*p/q).
nz1 = 0;
while ceil( ((Lx-1)*p+length(h)+nz1 )/q ) - delay < ceil(Lx*p/q)
    nz1 = nz1+1;
end
h = [h zeros(1,nz1)];

% ----  HERE'S THE CALL TO UPFIRDN  ----------------------------
y = upfirdn(x,h,p,q);

% Get rid of trailing and leading data so input and output signals line up
% temporally:
Ly = ceil(Lx*p/q);  % output length
% Ly = floor((Lx-1)*p/q+1);  <-- alternately, to prevent "running-off" the
%                                data (extrapolation)
if isvect
    y(1:delay) = [];
    y(Ly+1:end) = [];
else
    y(1:delay,:) = [];
    y(Ly+1:end,:) = [];
end

h([1:nz (end-nz1+1):end]) = [];  % get rid of leading and trailing zeros 
                                 % in case filter is output

%-------------------------------------------------------------------------
function [x, tx] = removeMissingTime(x,tx)
idx = isnan(tx);
if ~isempty(idx)
  tx(idx) = [];
  if isvector(x)
    x(idx) = [];
  else
    x(idx,:) = [];
  end
end

%-------------------------------------------------------------------------
function y = matInterp1(tin, xin, tout, method)

% by default specify tout as a column to obtain column matrix output
tout = tout(:);

if isvector(xin)
  if isrow(xin);
    % preserve orientation of input vector
    tout = tout.';
  end
  % interpolate, excluding NaN
  idx = find(~isnan(xin));  
  y = vecInterp1(tin(idx), xin(idx), tout, method);
else
  % initialize matrix output
  nRows = size(tout,1);
  nCols = size(xin,2);
  y = zeros(nRows,nCols);
  
  % loop through each column of input x
  for col=1:nCols
    % interpolate, excluding NaN
    idx = find(~isnan(xin(:,col)));  
    y(:,col) = vecInterp1(tin(idx), xin(idx,col), tout, method);
  end
end

%-------------------------------------------------------------------------
function y = vecInterp1(tin, xin, tout, method)

% check sample times for duplicate entries
iDup = find(diff(tin)==0);

% copy indices to point to the repeated locations in xin/tin
iRepeat = 1 + iDup;

while ~isempty(iDup)
  % find the number of successive equal sample times
  numEqual = find(diff(iDup)~=1,1,'first');
  if isempty(numEqual)
    numEqual = numel(iDup);
  end
  
  % replace leading x with mean value of all duplicates
  xSelect = xin(iDup(1) + (0:numEqual));
  xMean = mean(xSelect(~isnan(xSelect)));
  xin(iDup(1)) = xMean;

  % move to next block of conflicting sample times
  iDup = iDup(2+numEqual:end);
end

% remove duplicates
xin(iRepeat) = [];
tin(iRepeat) = [];

% call interp
y = interp1(tin, xin, tout, method, 'extrap');

%-------------------------------------------------------------------------
function [x, tx] = getSamples(x, tx)

validateattributes(x, {'numeric'},{'2d'}, ...
    'resample','X',1);
  
if isvector(x) && numel(x)~=numel(tx)
  error(message('signal:resample:TimeVectorMismatch'));
end

if ~isvector(x) && size(x,1)~=numel(tx);
  error(message('signal:resample:TimeRowMismatch'));
end

validateattributes(tx, {'numeric'},{'real','vector'}, ...
    'resample','Tx',2);

[x, tx] = removeMissingTime(x, tx);

% for efficiency, place samples in time-sorted order
[tx, idx] = sort(tx);
if isvector(x)
  % handle row vector input
  x = x(idx);
else
  x = x(idx,:);
end

% check finite behavior after sorting and NaN removal
validateattributes(tx, {'numeric'},{'finite'}, ...
    'resample','Tx',2);  

%-------------------------------------------------------------------------
function validateFs(fs)
validateattributes(fs, {'numeric'},{'real','finite','scalar','positive'}, ...
    'resample','Fs',3);

%-------------------------------------------------------------------------
function validateResampleRatio(p, q)
validateattributes(p, {'numeric'},{'integer','positive','finite','scalar'}, ...
    'resample','P');
validateattributes(q, {'numeric'},{'integer','positive','finite','scalar'}, ...
    'resample','Q');

%-------------------------------------------------------------------------
function [p, q] = getResampleRatio(t, fs)

% compute the average sampling interval
tsAvg = (t(end)-t(1))/(numel(t)-1);

% get a rational approximation of the ratio of the desired to the average
% sample rate
[p, q] = rat(tsAvg*fs,.01);

if p < 2
  % sample rate too small for input
  p = 1;
  q = round(1/(tsAvg*fs));
end

%-------------------------------------------------------------------------
function [method, arglist] = getInterpMethod(arglist)

method = '';
supportedMethods = {'linear','pchip','spline'};

iFound = 0;

for i=1:numel(arglist)
  if ischar(arglist{i})
    method = validatestring(arglist{i},supportedMethods,'resample','METHOD');
    iFound = i;
  end
end

if iFound
  arglist(iFound) = [];
end

