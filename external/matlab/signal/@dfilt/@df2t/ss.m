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

% Construct a df2 object
df2obj = dfilt.df2(Hd.Numerator,Hd.Denominator);

% Get its state-space representation
[a2,b2,c2,D] = ss(df2obj);

% Transpose a; b with c to get the df2t state-space representation
A = a2.';
% Force uniformity with empty state matrix.
if isempty(A)
  A = [];
  B = zeros(0,1);
  C = zeros(1,0);
else
  B = c2.';
  C = b2.';
end
