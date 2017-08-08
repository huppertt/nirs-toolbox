function [A,B,C,D] = ss(Hd)
%SS  Convert to state-space.
%   [A,B,C,D] = SS(Hd) converts discrete-time filter Hd to state-space
%   representation given by 
%     x(k+1) = A*x(k) + B*u(k)
%     y(k)   = C*x(k) + D*u(k)
%   where x is the state vector, u is the input vector, and y is the output
%   vector. 
%
%   See also DFILT.
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2002 The MathWorks, Inc.

A = Hd.A;
% Force uniformity with empty state matrix.
if isempty(A)
  A = [];
  B = zeros(0,1);
  C = zeros(1,0);
else
  B = Hd.B;
  C = Hd.C;
end
D = Hd.D;

% Do error checking on the consistency of the output.
error(abcdchk(A,B,C,D));