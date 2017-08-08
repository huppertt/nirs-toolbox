function [a,b] = statespaceab(Hd)
%STATESPACEAB A and B matrices of statespace realization.
%   [A, B] = statespaceab(Hd) returns the A and B matrices of the statespace
%   realization for auto-regressive lattice structures, which includes
%   Lattice AR, Lattice allpass, and Lattice ARMA.
  
%   Author(s): R. Losada, T. Bryan
%   Copyright 1988-2002 The MathWorks, Inc.

k = Hd.Lattice;

% Make sure k is a row
k = k(:).';
% Record the length of k
N = length(k);
% Generate a upper triangular matrix of ones
lowtriag = triu(ones(N),-1);
% Multiply each row except the first by conj(k(1:end-1))
for n = 1:N-1,
    lowtriag(n+1,:) = k(n)'.*lowtriag(n+1,:);
end
% Multiply each column by k
ktemp = k(ones(size(lowtriag,2),1),:);
lowtriag = lowtriag.*ktemp;
% Form the matrices a, b, c and d
a = [zeros(1,N);eye(N-1) zeros(N-1,1)] - lowtriag;
b = [1;k(1:end-1)'];
