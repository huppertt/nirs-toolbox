function P = classProbFromCost(C)
%
%   Example:
%      C = [0 2 1; 1 0 1; 1 2 0];
%      classProbFromCost(C)
%      % Compare with:
%      sum(C,2)/sum(C(:))
%
%   Copyright 2010 The MathWorks, Inc.


if isempty(C) || (~isnumeric(C) && ~islogical(C)) || ~ismatrix(C)
    error(message('stats:classreg:learning:classProbFromCost:BadCType'));
end

[K,M] = size(C);
if K~=M
    error(message('stats:classreg:learning:classProbFromCost:BadCSize'));
end

if any(diag(C)~=0)
    error(message('stats:classreg:learning:classProbFromCost:CNotDiagZero'));
end

if any(any(C<0))
    error(message('stats:classreg:learning:classProbFromCost:NegativeC'));
end

if isscalar(C)
    P = 1;
    return;
end

if any(all(C==0,2))
    error(message('stats:classreg:learning:classProbFromCost:CWithZeroRow'));
end

P = zeros(K,1);

% Nothing to scale for one class
if K==1
    P = 1;
end

% Unambiguous conversion for two classes
if K==2
    P(1) = C(1,2);
    P(2) = C(2,1);
    P = P/sum(P);
    return;
end

% Try to solve the linear system for multiple classes
N = (K-1)*K/2;
A = zeros(N,K);
Nblock = K-1;
Nfilled = 0;
for k=1:K-1
    A(Nfilled+1:Nfilled+Nblock,k)       = C(k+1:end,k);
    A(Nfilled+1:Nfilled+Nblock,k+1:end) = -diag(C(k,k+1:end));
    Nfilled = Nfilled + Nblock;
    Nblock = Nblock - 1;
end

% Get class probabilities for the solved system
if rank(A)<K
    % If there is a non-trivial solution, get it
    P = null(A);
else
    % Otherwise average the cost per class
    P = sum(C,2);
end

% Normalize
P = P/sum(P);

end
