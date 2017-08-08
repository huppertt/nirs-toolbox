function [z,p,H0,B,A] = ellipap2(N,Ap,As)
% Analog lowpass elliptic filter prototype.

%   Author(s): S. Orfanidis
%   Copyright 2006-2010 The MathWorks, Inc.


% Usage: [z,p,H0,B,A] = ellipap2(N,Ap,As)
%
% N  = filter order
% Ap = passband attenuation in dB
% As = stopband attenuation in dB
%
% z  = vector of normalized filter zeros (in units of the passband frequency Wp = 2*pi*fp)
% p  = vector of normalized filter poles
% H0 = DC gain factor
% B  = matrix whose rows are the first- and second-order numerator coefficients
% A  = matrix whose rows are the first- and second-order denominator coefficients
%
%        the gain factor g returned by ELLIPAP is related to the dc gain by
%        g = abs(H0*prod(p)/prod(z))
%
%        N = 2*L+r, r = mod(N,2), L = floor(N/2) = no. second-order sections
%
%        length(p) = N, length(z) = 2*L
%
%        transfer function: H(s) = H0 * [1/(1-s/p0)]^r * Prod_{i=1}^L [(1-s/zi)(1-s/zi^*)]/[(1-s/pi)(1-s/pi^*)]
%
%        normalized s-plane variable: s = j*Om/Om_p

Gp = 10^(-Ap/20);                                      % passband gain
ep = sqrt(10^(Ap/10) - 1);                             % ripple factors
es = sqrt(10^(As/10) - 1);

k1 = ep/es;                                 
k = ellipdeg(N,k1);                                    % solve degree equation

if k >= 1
  error(message('signal:ellipap2:FilterOrderTooLarge'));
end

L = floor(N/2); r = mod(N,2);                          % L is the number of second-order sections
i = (1:L)'; ui = (2*i-1)/N; zeta_i = cde(ui,k);        % zeros of elliptic rational function

j = complex(0,1);

z = j./(k*zeta_i);                                     % filter zeros = poles of elliptic rational function

v0 = -1i*asne(j/ep, k1)/N;                              % solution of sn(1i*v0*N*K1,k1) = j/ep

p = j*cde(ui-1i*v0, k);                                 % filter poles
p0 = j*sne(1i*v0, k);                                   % first-order pole, needed when N is odd

B = [ones(L,1), -2*real(1./z), abs(1./z).^2];          % second-order numerator sections
A = [ones(L,1), -2*real(1./p), abs(1./p).^2];          % second-order denominator sections

if r==0,                                               % prepend first-order sections
  B = [Gp, 0, 0; B];                                   % dc gain is Gp for N even
  A = [1, 0, 0; A];
else
  B = [1, 0, 0; B];                                    % dc gain is 1 for N odd
  A = [1, -real(1/p0), 0; A];
end

z = cplxpair([z; conj(z)]);                            % append conjugate zeros
p = cplxpair([p; conj(p)]);                            % append conjugate poles
if r==1, p = [p; p0]; end                              % append first-order pole when N is odd

H0 = Gp^(1-r);                                         % dc gain

