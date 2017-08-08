function [b,a] = prony(h, nb ,na)
%PRONY Prony's method for time-domain IIR filter design.
%   [B,A] = PRONY(H, NB, NA) finds a filter with numerator order
%   NB, denominator order NA, and having the impulse response in
%   vector H.   The IIR filter coefficients are returned in
%   length NB+1 and NA+1 row vectors B and A, ordered in
%   descending powers of Z.  H may be real or complex.
%
%   If the largest order specified is greater than the length of H,
%   H is padded with zeros.
%
%   % Example:
%   %   Fit an IIR model to an impulse response of a lowpass filter.
%
%   [b,a] = butter(4,0.2);
%   impulseResp = impz(b,a);                % obtain impulse response
%   denOrder=4; numOrder=4;                 % system function of order 4
%   [Num,Den]=prony(impulseResp,numOrder,denOrder);
%   subplot(211);                           % impulse response and input
%   stem(impz(Num,Den,length(impulseResp)));   
%   title('Impulse Response with Prony Design');
%   subplot(212);
%   stem(impulseResp); title('Input Impulse Response');
%
%   See also STMCB, LPC, BUTTER, CHEBY1, CHEBY2, ELLIP, INVFREQZ.

%   Author(s): L. Shure, 5-17-88
%              L. Shure, 12-17-90, revised
%   Copyright 1988-2013 The MathWorks, Inc.

%   References:
%     [1] T.W. Parks and C.S. Burrus, Digital Filter Design,
%         John Wiley and Sons, 1987, p226.

signal.internal.sigcheckfloattype(h,'','prony','H');
M = signal.internal.sigcasttofloat(nb,'double','prony','NB','allownumeric');
N = signal.internal.sigcasttofloat(na,'double','prony','NA','allownumeric');

K = length(h) - 1;

if K <= max(M,N)      % zero-pad input if necessary
    K = max(M,N)+1;
    h(K+1) = 0;
end
c = h(1);
if c==0    % avoid divide by zero
    c=1;
end
H = toeplitz(h/c,[1 zeros(1,K)]);
% K+1 by N+1
if (K > N)
    H(:,(N+2):(K+1)) = [];
end
% Partition H matrix
H1 = H(1:(M+1),:);	% M+1 by N+1
h1 = H((M+2):(K+1),1);	% K-M by 1
H2 = H((M+2):(K+1),2:(N+1));	% K-M by N
a = [1; -H2\h1].';
b = c*a*H1.';

