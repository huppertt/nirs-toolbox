function [R,sigma] = corrcov(C,nocheck)
%CORRCOV Compute correlation matrix from covariance matrix.
%   R = CORRCOV(C) computes the correlation matrix R that corresponds to the
%   covariance matrix C, by standardizing each row and column of C using the
%   square roots of the variances (diagonal elements) of C.  C is square,
%   symmetric, and positive semi-definite.  The correlation for a constant
%   variable (zero diagonal element of C) is undefined.
%
%   [R,SIGMA] = CORRCOV(C) computes the vector of standard deviations SIGMA
%   from the diagonal elements of C.
%
%   See also COV, CORR, CORRCOEF, CHOLCOV.

%   R = CORRCOV(C,1) computes the correlation matrix R without checking that C
%   is a valid covariance matrix.

%   Copyright 2007 The MathWorks, Inc. 


% Check square, symmetric, positive semidefinite.
if nargin < 2
    [T,p] = cholcov(C);
    if p ~= 0
        error(message('stats:corrcov:BadC'));
    end
end

[m,n] = size(C);
sigma = sqrt(diag(C)); % sqrt first to avoid under/overflow
R = bsxfun(@rdivide,C,sigma); R = bsxfun(@rdivide,R,sigma'); % R = C ./ sigma*sigma';

% Fix up possible round-off problems, while preserving NaN: put exact 1 on the
% diagonal, and limit off-diag to [-1,1]
t = find(abs(R) > 1); R(t) = R(t)./abs(R(t));
R(1:m+1:end) = sign(diag(R));
