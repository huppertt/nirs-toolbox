function [odata,b] = interp(idata,r,l,cutoff)
%INTERP Resample data at a higher rate using lowpass interpolation.
%   Y = INTERP(X,R) resamples the sequence in vector X at R times
%   the original sample rate.  The resulting resampled vector Y is
%   R times longer, LENGTH(Y) = R*LENGTH(X).
%
%   A symmetric filter, B, allows the original data to pass through
%   unchanged and interpolates between so that the mean square error
%   between them and their ideal values is minimized.
%   Y = INTERP(X,R,L,CUTOFF) allows specification of arguments
%   L and CUTOFF which otherwise default to 4 and .5 respectively.
%   2*L is the number of original sample values used to perform the
%   interpolation.  Ideally L should be less than or equal to 10.
%   The length of B is 2*L*R+1. The signal is assumed to be band
%   limited with cutoff frequency 0 < CUTOFF <= 1.0. 
%   [Y,B] = INTERP(X,R,L,CUTOFF) returns the coefficients of the
%   interpolation filter B.  
%
%   % Example:
%   %   Interpolate a signal by a factor of four.
%
%   t = 0:0.001:.029;                       % Time vector
%   x = sin(2*pi*30*t) + sin(2*pi*60*t);    % Original Signal
%   y = interp(x,4);                        % Interpolated Signal
%   subplot(211);
%   stem(x);
%   title('Original Signal');
%   subplot(212);
%   stem(y); 
%   title('Interpolated Signal');
%
%   See also DECIMATE, RESAMPLE, UPFIRDN.

%   Author(s): L. Shure, 5-14-87
%   	         L. Shure, 6-1-88, 12-15-88, revised
%   Copyright 1988-2014 The MathWorks, Inc.

%   References:
%     "Programs for Digital Signal Processing", IEEE Press
%     John Wiley & Sons, 1979, Chap. 8.1.

if nargin < 3
   l = 4;
end
if nargin < 4
   cutoff = .5;
end

% Check the input data type. Single precision is not supported.
try
    chkinputdatatype(idata,r,l,cutoff);
catch ME
    throwAsCaller(ME);
end

if l < 1 || r < 1 || cutoff <= 0 || cutoff > 1
   error(message('signal:interp:InvalidRange'))
end
if abs(r-fix(r)) > eps
   error(message('signal:interp:MustBeInteger'))
end
if 2*l+1 > length(idata)
	s = int2str(2*l+1);
	error(message('signal:interp:InvalidDimensions', s));
end

% ALL occurrences of sin()/() are using the sinc function for the
% autocorrelation for the input data. They should all be changed
% consistently if they are changed at all.

% calculate AP and AM matrices for inversion
s1 = toeplitz(0:l-1) + eps;
s2 = hankel(2*l-1:-1:l);
s2p = hankel([1:l-1 0]);
s2 = s2 + eps + s2p(l:-1:1,l:-1:1);
s1 = sin(cutoff*pi*s1)./(cutoff*pi*s1);
s2 = sin(cutoff*pi*s2)./(cutoff*pi*s2);
ap = s1 + s2;
am = s1 - s2;

% Compute matrix inverses using Cholesky decomposition for more robustness
U = chol(ap);
ap = inv(U)*inv(U).';
U = chol(am);
am = inv(U)*inv(U).';

% now calculate D based on INV(AM) and INV(AP)
d = zeros(2*l,l);
d(1:2:2*l-1,:) = ap + am;
d(2:2:2*l,:) = ap - am;

% set up arrays to calculate interpolating filter B
x = (0:r-1)/r;
y = zeros(2*l,1);
y(1:2:2*l-1) = (l:-1:1);
y(2:2:2*l) = (l-1:-1:0);
X = ones(2*l,1);
X(1:2:2*l-1) = -ones(l,1);
XX = eps + y*ones(1,r) + X*x;
y = X + y + eps;
h = .5*d'*(sin(pi*cutoff*XX)./(cutoff*pi*XX));
b = zeros(2*l*r+1,1);
b(1:l*r) = h';
b(l*r+1) = .5*d(:,l)'*(sin(pi*cutoff*y)./(pi*cutoff*y));
b(l*r+2:2*l*r+1) = b(l*r:-1:1);

% use the filter B to perform the interpolation
[m,n] = size(idata);
nn = max([m n]);
if nn == m
   odata = zeros(r*nn,1);
else
   odata = zeros(1,r*nn);
end
odata(1:r:nn*r) = idata;
% Filter a fabricated section of data first (match initial values and first derivatives by
% rotating the first data points by 180 degrees) to get guess of good initial conditions
% Filter length is 2*l*r+1 so need that many points; can't duplicate first point or
% guarantee a zero slope at beginning of sequence
od = zeros(2*l*r,1);
od(1:r:(2*l*r)) = 2*idata(1)-idata((2*l+1):-1:2);
[od,zi] = filter(b,1,od); %#ok
[odata,zf] = filter(b,1,odata,zi);
odata(1:(nn-l)*r) = odata(l*r+1:nn*r);

% make sure right hand points of data have been correctly interpolated and get rid of
% transients by again matching end values and derivatives of the original data
if nn == m
	od = zeros(2*l*r,1);
else
	od = zeros(1,2*l*r);
end
od(1:r:(2*l)*r) = 2*idata(nn)-(idata((nn-1):-1:(nn-2*l)));
od = filter(b,1,od,zf);
odata(nn*r-l*r+1:nn*r) = od(1:l*r);
