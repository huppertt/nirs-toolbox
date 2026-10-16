function dfe = satterthwaite_dfe(w, wXf)
% SATTERTHWAITE_DFE Fast Satterthwaite degrees of freedom estimate
%
%   dfe = nirs.math.satterthwaite_dfe(w, wXf)
%
%   Computes the Satterthwaite approximation to effective degrees of freedom
%   for robust regression without forming the full n×n hat matrix H.
%
%   This is equivalent to:
%       H = diag(w) - wXf*pinv(wXf'*wXf)*wXf';
%       HtH = H'*H;
%       dfe = trace(H'*H)^2 / trace((H'*H)^2);
%   but runs in O(n*k^2) instead of O(n^2) where k = size(wXf,2).
%
%   Inputs:
%       w    - [n×1] weights vector (e.g. from robust regression)
%       wXf  - [n×k] weighted, filtered design matrix (already multiplied by w)
%
%   Output:
%       dfe  - scalar effective degrees of freedom
%
%   The key insight: H = diag(w) - U*U' where U comes from SVD of wXf.
%   All required traces decompose into products of k×k matrices.
%
%   See also: nirs.math.ar_irls_fast, nirs.math.ar_irls

    w = w(:);  % ensure column

    % Economy SVD of weighted design matrix
    [U, ~, ~] = svd(wXf, 'econ');
    k = size(U, 2);

    % Leverage values: h_i = sum_j U_{ij}^2
    h = sum(U.^2, 2);

    % trace(H'H) = ||H||_F^2 where H = diag(w) - U*U'
    % = sum(w^2) - 2*sum(w.*h) + trace(U'U) = sum(w^2) - 2*sum(w.*h) + k
    trHtH = sum(w.^2) - 2 * sum(w .* h) + k;

    % trace((H'H)^2) = ||H^2||_F^2 via O(n*k^2) formula
    %
    % Since H is symmetric: H^2 = H'H, so we need trace(H^4).
    % H^2_{ii} = w_i^2 - (2*w_i - 1)*h_i
    % H^2_{ij} = -(w_i + w_j - 1)*P_{ij}  for i != j, where P = U*U'
    %
    % ||H^2||_F^2 = sum(diag_H2.^2) + sum_{i!=j} (w_i+w_j-1)^2 * P_{ij}^2
    %
    % The off-diagonal sum decomposes using k×k matrices:
    %   M_v = U' * diag(v) * U
    % and the identity sum_{ij} g(w_i)*h(w_j)*P_{ij}^2 = sum(M_g .* M_h, 'all')

    diag_H2 = w.^2 - (2*w - 1) .* h;

    % k×k matrices for trace decomposition
    Mw  = U' * bsxfun(@times, U, w);       % U'*diag(w)*U
    Mw2 = U' * bsxfun(@times, U, w.^2);    % U'*diag(w^2)*U
    Ik  = eye(k);                           % U'*U = I_k

    % (w_i + w_j - 1)^2 = w_i^2 + w_j^2 + 1 + 2*w_i*w_j - 2*w_i - 2*w_j
    % Full sum over ALL (i,j) including diagonal:
    offdiag_full = 2*trace(Mw2) ...               % w_i^2 + w_j^2 terms
                 + k ...                           % constant 1 term
                 + 2*sum(Mw(:).^2) ...             % 2*w_i*w_j term
                 - 4*trace(Mw);                    % -2*w_i - 2*w_j terms

    % Subtract diagonal contribution (i==j): (2*w_i - 1)^2 * h_i^2
    offdiag = offdiag_full - sum(((2*w - 1) .* h).^2);

    trH4 = sum(diag_H2.^2) + offdiag;

    dfe = trHtH^2 / max(trH4, eps);
end
