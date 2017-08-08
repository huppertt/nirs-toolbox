%TRIDISOLVE  Solve A*x = b where A is a square, symmetric tridiagonal matrix.
%   X = TRIDISOLVE(E,D,B,N) is the solution to the system
%     A*X = B
%   where A is an N-by-N symmetric, real, tridiagonal matrix given by
%     A = diag(E,-1) + diag(D,0) + diag(E,1)
%   Algorithm from Golub and Van Loan, "Matrix Computations", 2nd Edition, 
%   p.156.
%   Assumes A is non-singular.  If A is singular, a warning is issued and
%   results may be inaccurate.
%
%   Called by DPSS.
%
%   See also TRIDIEIG.

%   MEX-file implementation: T. Krauss, D. Orofino, 4/22/97
%   Copyright 1988-2002 The MathWorks, Inc.

%   MEX-file.

