function w = sne(u,k,tol)
%SNE  sn elliptic function with normalized complex argument.
%
% Usage: w = sne(u,k,tol)  (e.g., tol=1e-8)
%        w = sne(u,k,M)    (M=integer)
%        w = sne(u,k)      (equivalent to tol=eps)
%
% u = arbitrary vector of complex numbers on the u-plane 
% k = elliptic modulus (0 <= k < 1)
% tol = tolerance, e.g., tol=1e-8, default is tol = eps
%
% M = use a fixed number of Landen iterations, typically, M = 4-5
%
% w = the value of sn(u*K,k), w has the same size as u
%
% Notes: u is in units of the quarterperiod K. The conventional sn is computed by 
%        sn(u*K,k) = sne(u,k) ==> sn(v,k) = sne(v/K,k)
%
%        K = K(k), K' = K'(k) = K(k'), k' = sqrt(1-k^2)
%
%        k may not be 1, because K=Inf, but note that sn(u,1) = tanh(u)
%
%        it uses the Landen/Gauss transformation of ascending moduli
%        to build the answer from w = sin(u*pi/2)
%
%        the tolerance is that of computing the Landen vector v = landen(k,tol)
%
%        the ratio R=K'/K determines the pattern      
%        of zeros and poles of the sn function      (pole)N ---- D    u=j*R ---- u=1+j*R
%        within the SCDN fundamental rectangle,           |      |       |        |
%        the pole at corner N is u = j*R,                 |      |       |        |
%        the zero at corner S is u = 0              (zero)S ---- C      u=0 ---- u=1
%     
%        mappings around the S -> C -> D -> N path:
%             S -> C, 0<=t<=1, u = t        ==>    0 <= w <= 1     (passband)
%             C -> D, 0<=t<=1, u = 1+j*t*R  ==>    1 <= w <= 1/k   (transition)
%             D -> N, 0<=t<=1, u = 1-t+j*R  ==>  1/k <= w <= Inf   (stopband)
%
%        SNE and ASNE are inverses of each other.
%
%   See also LANDEN, CDE, ACDE, ELLIPK, ELLIPDEG, ELLIPJ, ELLIPKE.
        
%   Author(s): S. Orfanidis
%   Copyright 2009 The MathWorks, Inc.


if k==1, error(message('signal:sne:InvalidParam')); end
if nargin==2, tol=eps; end

v = landen(k,tol);                     % get Landen vector of descending moduli

w = sin(u*pi/2);

for n = length(v):-1:1,                % ascending Landen/Gauss transformation
   w = (1+v(n))*w ./ (1+v(n)*w.^2);
end


% [EOF]
