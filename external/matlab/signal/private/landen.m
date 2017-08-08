function v = landen(k,tol)
%LANDEN   Landen transformations of an elliptic modulus.
%
% Usage: v = landen(k,tol)     (e.g., tol = 1e-10)  
%        v = landen(k,M)       (M = integer)
%        v = landen(k)         (uses tol = eps)
%
% k  = starting elliptic modulus (must be 0 <= k <= 1)
% tol = tolerance, e.g., tol = 10^(-10), default value is tol = eps = 2.22e-16
%
% M = use fixed number of Landen iterations resulting in length(v) = M+1, typically, M = 4-5
%
% v  = Landen vector of descending-magnitude moduli starting with k
%
% Notes: The descending Landen/Gauss transformation is computed by the recurrence: 
%        v(n) = F(v(n-1)), for n = 2,3,..., 
%        initialized to v(1) = k, and iterated until v(n) < tol
%        where F(x) = [x/(1+sqrt(1-x^2))]^2
%  
%        the complete elliptic integral K(k) can be computed by K = prod(1 + v) * pi/2
% 
% examples: landen(0.95)      = [0.5241, 0.0801, 0.0016, 6.478e-7, 1.049e-13, 2.752e-27] (last k <= eps)
%           landen(0.95,1e-6) = [0.5241, 0.0801, 0.0016, 6.478e-7]                       (last k <= 1e-6)
%           landen(0.95,3)    = [0.5241, 0.0801, 0.0016]                                 (length set to 3)
%
%   See also ELLIPK, ELLIPDEG, ELLIPKE.

%   Author(s): S. Orfanidis
%   Copyright 2006 The MathWorks, Inc.

if nargin==1, tol=eps; end
if tol>=1, M=tol; end 

if k==0 || k==1, v=k; return; end  	% returns v=k, i.e., k=0 ==> v=0,  k=1 ==> v=1

v = [];

if tol<1,
   while k > tol,
      k = (k/(1+sqrt(1-k^2)))^2;
      v = [v; k];
   end
else
   for n=1:M, 
      k = (k/(1+sqrt(1-k^2)))^2;
      v = [v; k];
   end
end


% [EOF]
