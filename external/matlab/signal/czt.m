function g = czt(x, k, w, a)
%CZT  Chirp z-transform.
%   G = CZT(X, M, W, A) is the M-element z-transform of sequence X,
%   where M, W and A are scalars which specify the contour in the z-plane
%   on which the z-transform is computed.  M is the length of the transform,
%   W is the complex ratio between points on the contour, and A is the
%   complex starting point.  More explicitly, the contour in the z-plane
%   (a spiral or "chirp" contour) is described by
%       z = A * W.^(-(0:M-1))
%
%   The parameters M, W, and A are optional; their default values are 
%   M = length(X), W = exp(-1i*2*pi/M), and A = 1.  These defaults
%   cause CZT to return the z-transform of X at equally spaced points
%   around the unit circle, equivalent to FFT(X).
%
%   If X is a matrix, the chirp z-transform operation is applied to each
%   column.
%
%   % Example:
%   %   Use of czt to zoom in on a narrow-band section (100 to 150 Hz) of a 
%   %   filter's frequency response. 
%
%   %   Filter design
%   h = fir1(30,125/500,rectwin(31));   % Filter 
%
%   %   Establish frequency and CZT parameters
%   fs = 1000; f1 = 100; f2 = 150;      % In hertz 
%   m = 1024; 
%   w = exp(-j*2*pi*(f2-f1)/(m*fs)); 
%   a = exp(j*2*pi*f1/fs); 
%
%   %   Compute both the DFT and CZT of the filter
%   y = fft(h,1000); 
%   z = czt(h,m,w,a);
%
%   %   Create frequency vectors and compare the results
%   fy = (0:length(y)-1)'*1000/length(y);  
%   fz = ((0:m-1)'*(f2-f1)/m) + f1; 
%   plot(fy(1:500),abs(y(1:500))); axis([1 500 0 1.2]) 
%   title('FFT') 
%   figure 
%   plot(fz,abs(z)); axis([f1 f2 0 1.2]) 
%   title('CZT')
%
%   See also FFT, FREQZ.

%   Author(s): C. Denham, 1990.
%   	   J. McClellan, 7-25-90, revised
%   	   C. Denham, 8-15-90, revised
%   	   T. Krauss, 2-16-93, updated help
%   Copyright 1988-2013 The MathWorks, Inc.

%   References:
%     [1] Oppenheim, A.V. & R.W. Schafer, Discrete-Time Signal
%         Processing,  Prentice-Hall, pp. 623-628, 1989.
%     [2] Rabiner, L.R. and B. Gold, Theory and Application of
%         Digital Signal Processing, Prentice-Hall, Englewood
%         Cliffs, New Jersey, pp. 393-399, 1975.

[m, n] = size(x); oldm = m;
if m == 1, x = x(:); [m, n] = size(x); end

if nargin < 2
  k = length(x);
end

if nargin < 3
  % M cast to double to enforce precision rules
  Mtmp = signal.internal.sigcasttofloat(k,'double','czt','M','allownumeric');
  w = exp(-1i .* 2 .* pi ./ Mtmp); 
end
if nargin < 4
  a = 1;
end

% Checks if 'X', 'W' and 'A' are valid numeric data inputs
signal.internal.sigcheckfloattype(x,'','czt','X');
signal.internal.sigcheckfloattype(w,'','czt','W');
signal.internal.sigcheckfloattype(a,'','czt','A');

% Cast to enforce precision rules
k = signal.internal.sigcasttofloat(k,'double','czt','M','allownumeric');

if any([size(k) size(w) size(a)]~=1),
	error(message('signal:czt:InvalidDimensions'))
end

%------- Length for power-of-two fft.

nfft = 2^nextpow2(m+k-1);

%------- Premultiply data.

kk = ( (-m+1):max(k-1,m-1) ).';
kk2 = (kk .^ 2) ./ 2;
ww = w .^ (kk2);   % <----- Chirp filter is 1./ww
nn = (0:(m-1))';
aa = a .^ ( -nn );
aa = aa.*ww(m+nn);
y = x .* aa(:,ones(1,n));

%------- Fast convolution via FFT.

fy = fft(  y, nfft );
fv = fft( 1 ./ ww(1:(k-1+m)), nfft );   % <----- Chirp filter.
fy = fy .* fv(:,ones(1, n));
g  = ifft( fy );

%------- Final multiply.

g = g( m:(m+k-1), : ) .* ww( m:(m+k-1),ones(1, n) );

if oldm == 1, g = g.'; end

