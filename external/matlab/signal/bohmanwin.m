function w = bohmanwin(n)
%BOHMANWIN Bohman window.
%   BOHMANWIN(N) returns an N-point Bohman window in a column vector.
%
%   EXAMPLE:
%      N = 64; 
%      w = bohmanwin(N); 
%      plot(w); title('64-point Bohman window');
%
%   See also BARTHANNWIN, BARTLETT, BLACKMANHARRIS, FLATTOPWIN, 
%            NUTTALLWIN, PARZENWIN, RECTWIN, TRIANG, WINDOW.

%   Reference:
%     [1] fredric j. harris [sic], On the Use of Windows for Harmonic Analysis
%         with the Discrete Fourier Transform, Proceedings of the IEEE,
%         Vol. 66, No. 1, January 1978, Page 67, Equation 39.

%   Author(s): A. Dowd
%   Copyright 1988-2002 The MathWorks, Inc.

error(nargchk(1,1,nargin,'struct'));

[n,w,trivialwin] = check_order(n);
if trivialwin, return, end;

q = abs(linspace(-1,1,n));

% Forced end points to exactly zero
w = [ 0; ((1-q(2:end-1)).*cos(pi*q(2:end-1))+(1/pi)*sin(pi*q(2:end-1)))'; 0];


% [EOF]
