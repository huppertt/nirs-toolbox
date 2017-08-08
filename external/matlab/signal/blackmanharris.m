function w = blackmanharris(N,sflag)
%BLACKMANHARRIS Minimum 4-term Blackman-Harris window.
%   BLACKMANHARRIS(N) returns an N-point minimum 4-term Blackman-Harris 
%   window in a column vector.
%   BLACKMANHARRIS(N,SFLAG) generates the N-point Blackman-Harris window
%   using SFLAG window sampling. SFLAG may be either 'symmetric' or
%   'periodic'. By default, a symmetric window is returned.
%
%   EXAMPLE:
%      N = 32; 
%      w = blackmanharris(N); 
%      plot(w); title('32-point Blackman-Harris Window');
%
%   See also BARTHANNWIN, BARTLETT, BOHMANWIN, FLATTOPWIN, 
%            NUTTALLWIN, PARZENWIN, RECTWIN, TRIANG, WINDOW.

%   Reference:
%     [1] fredric j. harris [sic], On the Use of Windows for Harmonic 
%         Analysis with the Discrete Fourier Transform, Proceedings of 
%         the IEEE, Vol. 66, No. 1, January 1978

%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(1,2);    

if nargin < 2
    sflag = 'symmetric';
end

sflag = validatestring(sflag,{'symmetric','periodic'},'BLACKMANHARRIS','SFLAG');

% Cast to enforce Precision Rules
N = signal.internal.sigcasttofloat(N,'double','BLACKMANHARRIS','N',...
  'allownumeric');

% Check for valid window length (i.e., N < 0)
[N,w,trivialwin] = check_order(N);
if trivialwin, return, end;

% Coefficients obtained from page 65 of [1]
a = [0.35875 0.48829 0.14128 0.01168];
if strncmp(sflag,'p',1)
    x = (0:N-1)' * 2.0*pi/N;
else
    x = (0:N-1)'*2*pi/(N-1);
end
w = min4termwin(a,x);


% [EOF]
