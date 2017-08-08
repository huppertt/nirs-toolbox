function [z,p,k] = cheb2ap(n, rs)
%CHEB2AP Chebyshev Type II analog lowpass filter prototype.
%   [Z,P,K] = CHEB2AP(N,Rs) returns the zeros, poles, and gain
%   of an N-th order normalized analog prototype Chebyshev Type II
%   lowpass filter with Rs decibels of ripple in the stopband.
%   Chebyshev Type II filters are maximally flat in the passband.
%
%   % Example:
%   %   Design a 6-th order Chebyshev Type II analog low pass filter with
%   %   70dB of ripple in the stopband and display its frequency response.
%
%   [z,p,k]=cheb2ap(6,70);      % Lowpass filter prototype
%   [num,den]=zp2tf(z,p,k);     % Convert to transfer function form
%   freqs(num,den)              % Frequency response of analog filter  
%
%   See also CHEBY2, CHEB2ORD, BUTTAP, CHEB1AP, ELLIPAP.

%   Copyright 1988-2013 The MathWorks, Inc.

validateattributes(n,{'numeric'},{'scalar','integer','positive'},'cheb2ap','N');
validateattributes(rs,{'numeric'},{'scalar','nonnegative'},'cheb2ap','Rs');

% Cast to enforce precision rules
n = double(n);
rs = double(rs);

delta = 1/sqrt(10^(.1*rs)-1);
mu = asinh(1/delta)/n;
if (rem(n,2))
	m = n - 1;
	z = cos([1:2:n-2 n+2:2:2*n-1]*pi/(2*n))';
else
	m = n;
	z = cos((1:2:2*n-1)*pi/(2*n))';
end
z = (z - flipud(z))./2;
z = 1i./z;
% Organize zeros in complex pairs:
i = [1:m/2; m:-1:m/2+1];
z = z(i(:));

p = exp(1i*(pi*(1:2:2*n-1)/(2*n) + pi/2)).';
realp = real(p); realp = (realp + flipud(realp))./2;
imagp = imag(p); imagp = (imagp - flipud(imagp))./2;
p = complex(sinh(mu)*realp, cosh(mu)*imagp);
p = 1 ./ p;
k = real(prod(-p)/prod(-z));

