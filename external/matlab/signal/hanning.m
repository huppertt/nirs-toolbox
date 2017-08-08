function w = hanning(varargin)
%HANNING   Hanning window.
%   HANNING(N) returns the N-point symmetric Hanning window in a column
%   vector.  Note that the first and last zero-weighted window samples
%   are not included.
%
%   HANNING(N,'symmetric') returns the same result as HANNING(N).
%
%   HANNING(N,'periodic') returns the N-point periodic Hanning window,
%   and includes the first zero-weighted window sample.
%
%   NOTE: Use the HANN function to get a Hanning window which has the 
%          first and last zero-weighted samples. 
%
%   See also BARTLETT, BLACKMAN, BOXCAR, CHEBWIN, HAMMING, HANN, KAISER
%   and TRIANG.

%   Copyright 1988-2004 The MathWorks, Inc.

% Check number of inputs
error(nargchk(1,2,nargin,'struct'));

% Check for trivial order
[n,w,trivialwin] = check_order(varargin{1});
if trivialwin, return, end

% Select the sampling option
if nargin == 1,
   sflag = 'symmetric';
else
   sflag = lower(varargin{2});
end

% Allow partial strings for sampling options
allsflags = {'symmetric','periodic'};
sflagindex = find(strncmp(sflag, allsflags, length(sflag)));
if length(sflagindex)~=1         % catch 0 or 2 matches
   error(message('signal:hanning:InvalidEnum'));
end
sflag = allsflags{sflagindex};

% Evaluate the window
switch sflag,
case 'periodic'
   % Includes the first zero sample
   w = [0; sym_hanning(n-1)];
case 'symmetric'
   % Does not include the first and last zero sample
   w = sym_hanning(n);
end

%---------------------------------------------------------------------
function w = sym_hanning(n)
%SYM_HANNING   Symmetric Hanning window. 
%   SYM_HANNING Returns an exactly symmetric N point window by evaluating
%   the first half and then flipping the same samples over the other half.

if ~rem(n,2)
   % Even length window
   half = n/2;
   w = calc_hanning(half,n);
   w = [w; w(end:-1:1)];
else
   % Odd length window
   half = (n+1)/2;
   w = calc_hanning(half,n);
   w = [w; w(end-1:-1:1)];
end

%---------------------------------------------------------------------
function w = calc_hanning(m,n)
%CALC_HANNING   Calculates Hanning window samples.
%   CALC_HANNING Calculates and returns the first M points of an N point
%   Hanning window.

w = .5*(1 - cos(2*pi*(1:m)'/(n+1))); 

% [EOF] hanning.m
