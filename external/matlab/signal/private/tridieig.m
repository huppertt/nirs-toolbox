function x = tridieig(c,b,m1,m2,eps1);
%TRIDIEIG  Find a few eigenvalues of a tridiagonal matrix.
%   LAMBDA = TRIDIEIG(D,E,M1,M2).  D and E, two vectors of length N,
%   define a symmetric, tridiagonal matrix:
%      A = diag(E(2:N),-1) + diag(D,0) + diag(E(2:N),1)
%   E(1) is ignored.
%   TRIDIEIG(D,E,M1,M2) computes the eigenvalues of A with indices
%      M1 <= K <= M2.
%   TRIDIEIG(D,E,M1,M2,TOL) uses TOL as a tolerance.
%
%   See also TRIDISOLVE.

%   Author: C. Moler
%   Copyright 1988-2002 The MathWorks, Inc.

%   Mex-file version will be called if on the path.

%mbrealvector(c);
%mbrealvector(b);
%mbintscalar(m1);
%mbintscalar(m2);
%mbrealscalar(eps1)
if nargin < 5, eps1 = 0; end
n = length(c);
b(1) = 0;
beta = b.*b;
xmin = min(c(n) - abs(b(n)),min(c(1:n-1) - abs(b(1:n-1)) - abs(b(2:n))));
xmax = max(c(n) + abs(b(n)),max(c(1:n-1) + abs(b(1:n-1)) + abs(b(2:n))));
eps2 = eps*max(xmax,-xmin);
if eps1 <= 0, eps1 = eps2; end
eps2 = 0.5*eps1 + 7*eps2;

x0 = xmax;
x = zeros(n,1);
wu = zeros(n,1);
x(m1:m2) = xmax(ones(m2-m1+1,1));
wu(m1:m2) = xmin(ones(m2-m1+1,1));
z = 0;
for k = m2:-1:m1
   xu = xmin;
   for i = k:-1:m1
      if xu < wu(i)
         xu = wu(i);
         break
      end
   end
   if x0 > x(k), x0 = x(k); end
   while 1
      x1 = (xu + x0)/2;
      if x0 - xu <= 2*eps*(abs(xu)+abs(x0)) + eps1
         break
      end
      z = z + 1;
      a = 0;
      q = 1;
      for i = 1:n
         if q ~= 0
            s = beta(i)/q;
         else
            s = abs(b(i))/eps;
         end
         q = c(i) - x1 - s;
         a = a + (q < 0);
      end
      if a < k
         if a < m1
            xu = x1;
            wu(m1) = x1;
         else
            xu = x1;
            wu(a+1) = x1;
            if x(a) > x1, x(a) = x1; end
         end
      else
         x0 = x1;
      end
   end
   x(k) = (x0 + xu)/2;
end
x = x(m1:m2)';
