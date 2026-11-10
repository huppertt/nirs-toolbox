function r = sparseRank(X)
% Numerical rank via economy QR — works on both sparse and dense matrices
% without requiring full(X)
[~, R] = qr(X, 0);
d = abs(diag(R));
tol = max(size(X)) * eps(max(d));
r = sum(d > tol);
end
