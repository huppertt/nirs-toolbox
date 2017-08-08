function w = cde(u,k,tol)
%CDE  cd elliptic function with normalized complex argument.
%
% Usage: w = cde(u,k,tol)  (e.g., tol=1e-8)
%        w = cde(u,k,M)    (M=integer)
%        w = cde(u,k)      (equivalent to tol=eps)
%
% u = arbitrary vector of complex numbers on the u-plane
% k = elliptic modulus (0 <= k < 1)
% tol = tolerance, e.g., tol=1e-8, default is tol = eps
%
% M = use a fixed number of Landen iterations, typically, M = 4-5
%
% w = the value of cd(u*K,k), w has the same size as u
%
% Notes: u is in units of the quarterperiod K, thus, cd(z,k) = cde(z/K,k)
%
%        K = K(k), K' = K'(k) = K(k'), k' = sqrt(1-k^2)
%
%        k may not be 1, because K=Inf, but note that cd(u,1) = 1
%        
%        it uses the Landen/Gauss transformation of ascending moduli
%        to build the answer from w = cos(u*pi/2)
%
%        the tolerance is that of computing the Landen vector v = landen(k,tol)
%
%        the ratio R=K'/K determines the pattern      
%        of zeros and poles of the cd function      N ---- D(pole)    u=j*R ---- u=1+j*R
%        within the SCDN fundamental rectangle,     |      |             |        |
%        the pole at corner D is u = 1+j*R,         |      |             |        |
%        the zero at corner C is u = 1              S ---- C(zero)      u=0 ---- u=1
%     
%        mappings around the C -> S -> N -> D path:
%             C -> S, 0<=t<=1, u = 1-t    ==>    0 <= w <= 1     (passband)
%             S -> N, 0<=t<=1, u = j*t*R  ==>    1 <= w <= 1/k   (transition)
%             N -> D, 0<=t<=1, u = t+j*R  ==>  1/k <= w <= Inf   (stopband)
%
%        CDE and ACDE are inverses of each other.
%
%   See also LANDEN, SNE, ASNE, ELLIPK, ELLIPDEG, ELLIPJ, ELLIPKE.
        
%   Author(s): S. Orfanidis
%   Copyright 2006-2009 The MathWorks, Inc.

validateattributes(k,{'numeric'},{'scalar','real','nonempty','<',1,'>=',0},...
    'cde','k');

if nargin==2, tol=eps; end

v = landen(k,tol);                     % get Landen vector of descending moduli

w = cos(u*pi/2);

for n = length(v):-1:1,                % ascending Landen/Gauss transformation
   w = (1+v(n))*w ./ (1+v(n)*w.^2);
end


% [EOF]
