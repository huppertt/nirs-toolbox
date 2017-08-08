function [z,p,k] = buttap(n)
%BUTTAP Butterworth analog lowpass filter prototype.
%   [Z,P,K] = BUTTAP(N) returns the zeros, poles, and gain
%   for an N-th order normalized prototype Butterworth analog
%   lowpass filter.  The resulting filter has N poles around
%   the unit circle in the left half plane, and no zeros.
%
%   % Example:
%   %   Design a 9th order Butterworth analog lowpass filter and display
%   %   its frequency response.
%
%   [z,p,k]=buttap(9);          % Butterworth filter prototype
%   [num,den]=zp2tf(z,p,k);     % Convert to transfer function form
%   freqs(num,den)              % Frequency response of analog filter          
%
%   See also BUTTER, CHEB1AP, CHEB2AP, ELLIPAP.

%   Author(s): J.N. Little and J.O. Smith, 1-14-87
%   	   L. Shure, 1-13-88, revised
%   Copyright 1988-2002 The MathWorks, Inc.

validateattributes(n,{'numeric'},{'scalar','integer','positive'},'buttap','N');
% Cast to enforce precision rules
n = double(n);
% Poles are on the unit circle in the left-half plane.
z = [];
p = exp(1i*(pi*(1:2:n-1)/(2*n) + pi/2));
p = [p; conj(p)];
p = p(:);
if rem(n,2)==1   % n is odd
    p = [p; -1];
end
k = real(prod(-p));

