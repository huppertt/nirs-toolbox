function LL = fitlme_loglik(X, Y, Z, theta)
% Compute profile log-likelihood for a linear mixed-effects model at given theta
% LL = nirs.math.fitlme_loglik(X, Y, Z, theta)
%
% Lightweight version of the solveLME inner loop — returns only the
% profile log-likelihood.  Used for joint theta optimization across
% channels in the Kronecker decomposition of MixedEffects.
%
% Inputs:
%   X     - Fixed effects design matrix (nT x nX)
%   Y     - Response vector (nT x 1)
%   Z     - Random effects design matrix (nT x nZ)
%   theta - Random effects variance parameter (scalar)
%
% Output:
%   LL    - Profile log-likelihood value

if nargin < 4, theta = 0; end

nT = size(Y, 1);
nZ = size(Z, 2);

% Handle NaN rows
bad = any(isnan(X),2) | any(isnan(Y),2) | any(isnan(Z),2);
if any(bad)
    X(bad,:) = [];
    Y(bad,:) = [];
    Z(bad,:) = [];
    nT = size(Y, 1);
end

if nT == 0
    LL = -Inf;
    return;
end

% No random effects: simple OLS log-likelihood
if isempty(Z) || nZ == 0
    [Q, R] = qr(X, 0);
    beta = R \ (Q' * Y);
    resid = Y - X * beta;
    r2 = sum(resid.^2);
    if r2 <= 0, LL = -Inf; return; end
    LL = (-nT/2) * (1 + log(2*pi*r2/nT));
    return;
end

% Build Lambda and factor the random effects precision
Lambda = sqrt(exp(theta)) * speye(nZ);
Iq = speye(nZ);
sZ = sparse(Z);

a = Lambda' * (sZ'*sZ) * Lambda + Iq;
[R, p, S] = chol(a);
if p ~= 0
    LL = -Inf;
    return;
end

Q1 = ((X'*sZ*Lambda)*S) / R;
R1R1t = X'*X - Q1*Q1';

% Safe Cholesky with regularization
delta = eps(class(R1R1t));
if issparse(R1R1t)
    I = speye(size(R1R1t));
else
    I = eye(size(R1R1t));
end
p = 1; iter = 0;
while p ~= 0
    try
        [R1, p] = chol(R1R1t + delta*I, 'lower');
    catch
        p = 1;
    end
    delta = 2 * delta;
    iter = iter + 1;
    if iter > 1000, LL = -Inf; return; end
end

% Parameter estimates (only needed to compute residuals)
cDeltab = R' \ (S' * (Lambda' * Z' * Y));
cbeta = R1 \ (X' * Y - Q1 * cDeltab);
beta = R1' \ cbeta;
Deltab = S * (R \ (cDeltab - Q1' * beta));

% Residuals and error
resid = Y - X * beta - Z * (Lambda * Deltab);
r2 = sum(Deltab.^2) + sum(resid.^2);

if r2 <= 0, LL = -Inf; return; end

% Profile log-likelihood
LL = (-nT/2) * (1 + log(2*pi*r2/nT)) - sum(log(diag(R)));

end
