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

% Normalize the transfer function and get order of num and den
num = num./den(1); den = den./den(1);
N = length(num)-1; M = length(den) - 1;

if M == 0 & N == 0
    % Scalar case,
    A = [];
    B = zeros(0,1);
    C = zeros(1,0);
    D = num;
    return
end

if M==0,
    % FIR
    A = sysmatrix(den,N,M);
    B = [1;
        zeros(N-1,1)];
elseif N==0,
    % All pole
    A = [-den(2:end);
        eye(M-1) zeros(M-1,1)];
    B = [num(1);
        zeros(M-1,1)];
else
    % IIR
    % Build the system matrix A
    A = [zeros(1,N+M);
        eye(N-1) zeros(N-1,M+1);
        num(2:end) -den(2:end);
        zeros(M-1,N) eye(M-1) zeros(M-1,1)];
    
    % Build the input matrix B
    B = [1;
        zeros(N-1,1);
        num(1);
        zeros(M-1,1)];
end

% Build output matrix C
C = [num(2:end) -den(2:end)];

% Build direct feedthrough matrix D
D = num(1);

if isempty(A)
  A = [];
  B = zeros(0,1);
  C = zeros(1,0);
end
