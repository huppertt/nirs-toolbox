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

num = Hd.Numerator;
den = Hd.Denominator;

% Get order of num and den
N = length(num)-1; M = length(den) - 1;

% Normalize the transfer function
num = num./den(1); den = den./den(1);

if M <= 0 & N <= 0
    % Empty or Scalar case,
    A = []; 
    B = zeros(0,1);
    C = zeros(1,0);
    D = num;
    return
end

% Make sure num is at least as long as den
num = [num zeros(1,M-N)];
% Recalculate the num order
N = length(num)-1;

% Build the system matrix A
A = sysmatrix(den,N,M);
	
% Force uniformity with empty state matrix.
if isempty(A)
  A = [];
  B = zeros(0,1);
  C = zeros(1,0);
else
  % Build the input matrix B
  B = [1;
       zeros(N-1,1)];
	
  % Build output matrix C
  ct1 = [-num(1).*den(2:end)];
  ct2 = [num(2:end)];
  ct1 = [ct1 zeros(1,length(ct2)-length(ct1))];
  ct2 = [ct2 zeros(1,length(ct1)-length(ct2))];
  C = ct1+ct2;
end

% Build direct feedthrough matrix D
D = num(1);

