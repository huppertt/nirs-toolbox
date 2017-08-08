function k = ellipdeg(N,k1,tol)
%ELLIPDEG Solves the degree equation in analog elliptic filter design.
%
% Usage: k = ellipdeg(N,k1,tol);  
%        k = ellipdeg(N,k1,M);    (M = fixed number of Landen iterations)
%        k = ellipdeg(N,k1)       (uses tol=eps)
%
% N = analog filter order
% k1 = elliptic modulus for stopband band, that is, k1 = ep/ep1
% tol = tolerance, e.g., tol=1e-10, default is tol = eps
%
% M = uses a fixed number of Landen iterations, typically, M = 4-5
%
% k = elliptic modulus for transition band, that is, k = WB/W1
%
% Notes: solves the degree equation N*K'/K = K1'/K1 for k in terms of N,k1
%        it uses Jacobi's exact solution k' = (k1')^N * (Prod_i(sn(ui*K1',k1')))^4
%        when k1 is very small, it uses the function ELLIPDEG2 to avoid numerical
%        problems arising from the computation of sqrt(1-k1^2)
%        to solve for k1, given N,k, use the function ELLIPDEG1

%   Author(s): S. Orfanidis
%   Copyright 2006 The MathWorks, Inc.

if nargin==2, tol=eps; end

L = floor(N/2); 

i = 1:L;  ui = (2*i-1)/N;

kmin = 1e-6;

if k1 < kmin
   k = ellipdeg2(1/N,k1,tol);              % use nome method when k1 is too small
else                                       
   kc = sqrt(1-k1^2);			   % complement of k1 
   kp = kc^N * (prod(sne(ui,kc)))^4;       % complement of k
   k = sqrt(1-kp^2); 
end

% when k1 is near 1, then so is k, because k1 < k < 1, and if the calculated kp 

% becomes less than about eps, then k = sqrt(1-kp^2) would be inaccurate. 

% To avoid this, N may not be too large. The maximum usable value of N consistent

% with Matlab's floating point accuracy is when k becomes equal to about k=1-eps,

% which gives Nmax = log(q1)/log(q), with q1 = nome(k1) and q = nome(1-eps), or,

% approximately, Nmax = -3.86 * log(q1) 

%%
function k1 = ellipdeg2(N,k,tol,M)
%ellipdeg2 - solves the degree equation in analog elliptic filter design
%
% Usage: k1 = ellipdeg2(N,k,tol,M);
%        k1 = ellipdeg2(N,k,tol);    (use M=7 expansion terms)
%        k1 = ellipdeg2(N,k)         (uses tol=eps, M=7)
%
%        k = ellipdeg2(1/N, k1, tol)  implements the inverse operation
%
% N = analog filter order
% k = elliptic modulus for transition band
% tol = tolerance, e.g., tol=1e-10, default is tol = eps
% M = number of series expansion terms, default M=7
%
% k1 = elliptic modulus for stopband
%
% Notes: the degree equation N*K(k1)/K(k) = K'(k1)/K'(k) is solved for k1, given N,k
%        or, for k, given N,k1
%
%        It computes the nome q of k, then the nome of k1 from q1 = q^N, and calculates
%        k1 from the approximation k1 = 4*sqrt(q1) * (higher orders in q1)
%
%        The inverse of k1 = ellipdeg(N,k) can be calculated from k = ellipdeg(1/N,k1)
%        alternatively, k' = ellipdeg(N,k1'), where k1'=sqrt(1-k1^2), k' = sqrt(1-k^2)
%        
%        it calls ELLIPK to calculate the nome q of the modulus k

if nargin<=3, M=7; end
if nargin==2, tol=eps; end

[K,Kprime] = ellipk(k,tol);  q = exp(-pi*Kprime/K);

q1 = q^N;

m = 1:M;

k1 = 4 * sqrt(q1) * ((1 + sum(q1.^(m.*(m+1))))/(1 + 2*sum(q1.^(m.*m))))^2;


% [EOF]
