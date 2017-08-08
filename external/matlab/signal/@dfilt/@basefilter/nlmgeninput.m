function [vp,Vp] = nlmgeninput(this,M,L)
%NLMGENINPUT   Generate input signal for noisepsd and mfreqresp.
%   Author(s): R. Losada
%   Copyright 1988-2006 The MathWorks, Inc.

if nargin < 3, L = 1; end

if rem(M,2) == 0,
    phi = 2*pi*rand(M/2-1,L);           % Generation of the periodic
    phi = [zeros(1,L); phi; zeros(1,L); -phi(end:-1:1,:)];  % input signal.
else
    phi = 2*pi*rand((M-1)/2,L);           % Generation of the periodic
    phi = [zeros(1,L); phi;-phi(end:-1:1,:)];  % input signal.
end
Vp = exp(i*phi);
vp = ifft(Vp,'symmetric');

% [EOF]
