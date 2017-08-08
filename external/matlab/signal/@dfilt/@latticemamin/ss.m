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

k = Hd.Lattice;

% Make sure k is A column
k = k(:);

% Record the length of k
N = length(k);

% Handle the scalar case
if N==1 & k == 0,
    A = [];
    B = zeros(0,1);
    C = zeros(1,0);
    D = 1;
    return
end


% Generate A lower triangular matrix of ones
lowtriag = tril(ones(N-2));

% Multiply each row by conj(k(2:end-1))
ktemp = conj(k(2:end-1,ones(size(lowtriag,2),1)));
lowtriag = lowtriag.*ktemp;

% Multiply each column by k(1:end-2)
ktemp = k(1:end-2).'; % First use A row vector of k
ktemp = ktemp(ones(size(lowtriag,2),1),:);
if ~isempty(lowtriag), % Necessary if length(k) == 1
    lowtriag = lowtriag.*ktemp;
end

% Add A row and A column to the lower triang matrix plus an identity
lowtriag = eye(N-1) + [zeros(1,N-1);lowtriag zeros(N-2,1)];

% Form the system matrix A
A = [zeros(1,N);
    lowtriag zeros(N-1,1)];
    
% Force uniformity with empty state matrix.
if isempty(A)
  A = [];
  B = zeros(0,1);
  C = zeros(1,0);
else
  B = [1;conj(k(1:end-1))];
  C = k.';
end

D = 1;

