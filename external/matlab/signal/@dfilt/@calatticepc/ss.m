function [A,B,C,D] = ss(Hd)
%SS  Discrete-time filter to state-space conversion.
%   [A,B,C,D] = SS(Hd) converts discrete-time filter Hd to state-space
%   representation given by 
%     x(k+1) = A*x(k) + B*u(k)
%     y(k)   = C*x(k) + D*u(k)
%   where x is the state vector, u is the input vector, and y is the output
%   vector. 
%
%   See also DFILT.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

[A1,B1,C1,D1] = ss(Hd.privAllpass1);
[A2,B2,C2,D2] = ss(Hd.privAllpass2);

beta = Hd.beta;

A = [A1, zeros(size(A1,1),size(A2,2));
     zeros(size(A2,1),size(A1,2)), A2];

if isempty(A)
  A = [];
  B = zeros(0,1);
  C = zeros(1,0);
else
  B = (i/2)*[B1;B2];
  C = [conj(beta)*C1, -beta*C2];
end

D = (i/2)*(conj(beta)*D1 - beta*D2);
