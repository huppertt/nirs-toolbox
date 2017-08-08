function [z,p,k] = cheb1ap(n, rp)
%CHEB1AP Chebyshev Type I analog lowpass filter prototype.
%   [Z,P,K] = CHEB1AP(N,Rp) returns the zeros, poles, and gain
%   of an N-th order normalized analog prototype Chebyshev Type I
%   lowpass filter with Rp decibels of ripple in the passband.
%   Chebyshev Type I filters are maximally flat in the stopband.
%
%   % Example:
%   %   Design a 6-th order Chebyshev Type I analog low pass filter with 
%   %   3dB of ripple in the passbandand and display its frequency 
%   %   response.
%
%   [z,p,k]=cheb1ap(6,3);       % Lowpass filter prototype
%   [num,den]=zp2tf(z,p,k);     % Convert to transfer function form
%   freqs(num,den)              % Frequency response of analog filter  
%  
%   See also CHEBY1, CHEB1ORD, BUTTAP, CHEB2AP, ELLIPAP.

%   Copyright 1988-2013 The MathWorks, Inc.

validateattributes(n,{'numeric'},{'scalar','integer','positive'},'cheb1ap','N');
validateattributes(rp,{'numeric'},{'scalar','nonnegative'},'cheb1ap','Rp');

% Cast to enforce precision rules
n = double(n);
rp = double(rp);

epsilon = sqrt(10^(.1*rp)-1);
mu = asinh(1/epsilon)/n;
p = exp(1i*(pi*(1:2:2*n-1)/(2*n) + pi/2)).';
realp = real(p); realp = (realp + flipud(realp))./2;
imagp = imag(p); imagp = (imagp - flipud(imagp))./2;
p = complex(sinh(mu).*realp , cosh(mu).*imagp);
z = [];
k = real(prod(-p));
if ~rem(n,2)	% n is even so patch k
	k = k/sqrt((1 + epsilon^2));
end
