function w = barthannwin(N)
%BARTHANNWIN Modified Bartlett-Hanning window.
%   BARTHANNWIN(N) returns an N-point Modified Bartlett-Hanning window
%   in a column vector.
%
%   EXAMPLE:
%      N = 64;
%      w = barthannwin(N);
%      plot(w); title('64-point Modified Bartlett-Hanning window');
%
%   See also BARTLETT, BLACKMANHARRIS, BOHMANWIN, FLATTOPWIN,
%            NUTTALLWIN, PARZENWIN, RECTWIN, TRIANG, WINDOW.

%   Reference:
%     [1] Yeong Ho Ha and John A. Pearce, A New Window and Comparison
%         to Standard Windows, IEEE Transactions on Acoustics, Speech,
%         and Signal Processing, Vol. 37, No. 2, February 1999

%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(1,1);  

% Cast to enforce Precision Rules
N = signal.internal.sigcasttofloat(N,'double','barthannwin','N',...
  'allownumeric');

% Check for valid window length (i.e., N < 0)
[N,w,trivialwin] = check_order(N);
if trivialwin, return, end;

% Equation 6 from [1]
t = ((0:(N-1))/(N-1)-.5)'; % -.5 <= t <= .5
w = .62-.48*abs(t) + .38*cos(2*pi*t);


% [EOF]
