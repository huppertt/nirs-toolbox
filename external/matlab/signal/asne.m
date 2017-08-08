function u = asne(w,k,tol)
%ASNE  Inverse of sn elliptic function.
%
% Usage: u = asne(w,k,tol)  (e.g., tol=1e-8)
%        u = asne(w,k,M)    (M=integer)
%        u = asne(w,k)      (equivalent to tol=eps)
%
% w = arbitrary vector of complex numbers on the u-plane
% k = elliptic modulus (0 <= k < 1)
% tol = tolerance, e.g., tol=1e-8, default is tol = eps
% M = use a fixed number of Landen iterations, typically, M = 4-5
%
% u = the solution of sn(u*K,k) = w, u has the same size as w, u=0 if w=0
%
% Notes: u is in units of the quarterperiod K, so that w = sn(u*K,k) = sne(u,k) is inverted
%        by  u = asne(w,k)
%
%        K = K(k), K' = K'(k) = K(k'), k' = sqrt(1-k^2)
%                                                                                   2 | 1
%        u is reduced into the rectangle, -1<Re(u)<1, -R<Im(u)<R, quadrant mapping  -----, centered at u=0
%                                                                                   3 | 4
%        k may not be 1, because K=Inf, but note that sn(u,1) = tanh(u)
%        
%        it uses the property: sn(u*K,k) = cd(K*(1-u),k) ==> u = 1 - acde(u,k,tol)  
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

if k==1, error(message('signal:asne:InvalidParam')); end
if nargin==2, tol=eps; end

u = 1-acde(w,k,tol);    


% [EOF]
