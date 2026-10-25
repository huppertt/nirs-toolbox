function covb = build_covb_tensor(Xfall, C)
% BUILD_COVB_TENSOR Fast 4D cross-channel covariance tensor
%
%   covb = nirs.math.build_covb_tensor(Xfall, C)
%
%   Builds the 4D covariance tensor covb(:,:,i,j) from weighted filtered
%   design matrices (Xfall) and cross-channel covariance scalars (C).
%
%   This replaces the common pattern:
%       for i=1:nChan
%           for j=1:nChan
%               covb(:,:,i,j) = pinv(Xfall{i}'*Xfall{j}) * C(i,j);
%           end
%       end
%
%   Optimizations:
%     1. Batch Gram matrix via single BLAS call (no-NaN fast path)
%     2. Upper triangle only — fills both (i,j) and (j,i) per pair
%     3. Backslash instead of pinv for small well-conditioned matrices
%
%   Inputs:
%       Xfall - {nChan x 1} cell array, each nTime x nCond (may contain NaN)
%       C     - [nChan x nChan] cross-channel covariance scalar matrix
%
%   Output:
%       covb  - [nCond x nCond x nChan x nChan] covariance tensor
%
%   See also: nirs.math.ar_irls_fast, nirs.math.satterthwaite_dfe

    nChan = numel(Xfall);
    nCond = size(Xfall{1}, 2);
    nTime = size(Xfall{1}, 1);

    % Check if any Xfall has NaN (channels with missing timepoints)
    hasNaN = false;
    for i = 1:nChan
        if any(isnan(Xfall{i}(:)))
            hasNaN = true;
            break;
        end
    end

    covb = zeros(nCond, nCond, nChan, nChan);
    Ik = eye(nCond);

    if ~hasNaN
        % FAST PATH: batch Gram matrix via single BLAS call
        % Stack all weighted filtered design matrices into one matrix
        Wmat = zeros(nTime, nCond * nChan);
        for i = 1:nChan
            cols = (i-1)*nCond + (1:nCond);
            Wmat(:, cols) = Xfall{i};
        end

        % One matrix multiply replaces nChan^2 individual Xi'*Xj products
        GG = Wmat' * Wmat;  % (nCond*nChan) x (nCond*nChan)

        % Upper triangle only, backslash with pinv fallback
        for i = 1:nChan
            idx_i = (i-1)*nCond + (1:nCond);
            for j = i:nChan
                idx_j = (j-1)*nCond + (1:nCond);
                Gij = GG(idx_i, idx_j);
                Gji = GG(idx_j, idx_i);
                rhs_ij = Ik * C(i, j);
                rhs_ji = Ik * C(j, i);
                w = warning('off', 'MATLAB:singularMatrix');
                covb(:, :, i, j) = Gij \ rhs_ij;
                covb(:, :, j, i) = Gji \ rhs_ji;
                warning(w);
                if any(~isfinite(covb(:, :, i, j)), 'all')
                    covb(:, :, i, j) = pinv(Gij) * rhs_ij;
                end
                if any(~isfinite(covb(:, :, j, i)), 'all')
                    covb(:, :, j, i) = pinv(Gji) * rhs_ji;
                end
            end
        end
    else
        % SLOW PATH: per-pair NaN handling, upper triangle only
        for i = 1:nChan
            Xi = Xfall{i};
            for j = i:nChan
                Xj = Xfall{j};
                lstV = ~isnan(sum(Xi, 2) + sum(Xj, 2));
                Gij = Xi(lstV, :)' * Xj(lstV, :);
                Gji = Xj(lstV, :)' * Xi(lstV, :);
                rhs_ij = Ik * C(i, j);
                rhs_ji = Ik * C(j, i);
                w = warning('off', 'MATLAB:singularMatrix');
                covb(:, :, i, j) = Gij \ rhs_ij;
                covb(:, :, j, i) = Gji \ rhs_ji;
                warning(w);
                if any(~isfinite(covb(:, :, i, j)), 'all')
                    covb(:, :, i, j) = pinv(Gij) * rhs_ij;
                end
                if any(~isfinite(covb(:, :, j, i)), 'all')
                    covb(:, :, j, i) = pinv(Gji) * rhs_ji;
                end
            end
        end
    end
end
