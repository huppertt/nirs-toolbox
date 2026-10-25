function [stats, resid] = ar_irls_fast(d, X, Pmax, tune, nosearch, useGPU, singlePrecision, maxBICSearch)
% AR_IRLS_FAST - Parallelized AR(P)-IRLS regression
%
% Drop-in replacement for nirs.math.ar_irls with these optimizations:
%   1. parfor across channels (falls back to for-loop if no pool)
%   2. Satterthwaite DOF without forming n×n H matrix (O(n*k^2) vs O(n^2))
%   3. Vectorized cross-channel covariance
%   4. Capped BIC search + skip BIC after first IRLS iteration per channel
%
% Usage:
%   [stats, resid] = nirs.math.ar_irls_fast(d, X, Pmax, tune, nosearch, useGPU, singlePrecision, maxBICSearch)
%
%   maxBICSearch - (optional) Maximum AR orders to evaluate in BIC search.
%       0 = full search up to Pmax every iteration (default, same as ar_irls)
%       N>0 = search orders 1..N on first IRLS iteration, then direct fit
%             at found order on subsequent iterations (much faster)
%       Typical value: 3 (sufficient for most fNIRS data)
%
% See also: nirs.math.ar_irls

    warning('off', 'stats:statrobustfit:IterationLimit');

    if nargin < 4 || isempty(tune), tune = 4.685; end
    if nargin < 5, nosearch = false; end
    if nargin < 6, useGPU = false; end
    if nargin < 7, singlePrecision = false; end
    if nargin < 8 || isempty(maxBICSearch), maxBICSearch = 0; end

    nCond = size(X, 2);
    nChan = size(d, 2);
    nTime = size(d, 1);

    % Preallocate outputs as cell arrays for parfor compatibility
    beta_all = zeros(nCond, nChan);
    tstat_all = zeros(nCond, nChan);
    pval_all = zeros(nCond, nChan);
    ppos_all = zeros(nCond, nChan);
    pneg_all = zeros(nCond, nChan);
    P_all = zeros(nChan, 1);
    w_all = zeros(nTime, nChan);
    dfe_all = zeros(nChan, 1);
    covb_per = zeros(nCond, nCond, nChan);
    sigma2_all = zeros(nChan, 1);
    a_all = cell(nChan, 1);
    filter_all = cell(nChan, 1);
    R2_all = zeros(nChan, 1);
    logLik_all = zeros(nChan, 1);
    resid = nan(nTime, nChan);
    Xfall = cell(nChan, 1);

    % Precompute initial OLS for warm-start (shared across channels, small cost)
    % Note: X_pinv(:,lstValid)*y(lstValid) is equivalent to pinv(X(lstValid,:))*y(lstValid)
    % only when lstValid is all true (no NaN). For safety, we do pinv inside loop when NaNs exist.
    X_pinv = pinv(X);
    hasNoNans = ~any(isnan(d(:)));

    % Handle per-channel Pmax
    if length(Pmax) == 1
        Pmax_vec = repmat(Pmax, nChan, 1);
    else
        Pmax_vec = Pmax(:);
    end

    % ============================================================
    % MAIN CHANNEL LOOP - parallelized
    % ============================================================
    parfor i = 1:nChan
        y = d(:, i);

        if ~any(y)
            % Zero channel - fill with NaN
            Xfall{i} = nan(nTime, nCond);
            beta_all(:, i) = NaN;
            covb_per(:, :, i) = NaN;
            w_all(:, i) = NaN;
            a_all{i} = [];
            sigma2_all(i) = NaN;
            tstat_all(:, i) = NaN;
            pval_all(:, i) = NaN;
            ppos_all(:, i) = NaN;
            pneg_all(:, i) = NaN;
            filter_all{i} = [];
            R2_all(i) = NaN;
            logLik_all(i) = NaN;
            resid(:, i) = NaN;
            continue;
        end

        lstValid = ~isnan(y);
        lstInvalid = ~lstValid;
        nValid = sum(lstValid);

        % Initial OLS fit
        if hasNoNans
            B = X_pinv * y;
        else
            B = pinv(X(lstValid,:)) * y(lstValid);
        end
        B0 = 1e6 * ones(size(B));

        p_i = Pmax_vec(i);
        iter = 0;
        maxiter = 10;
        f = 1;
        Xf = X;
        yf = y;
        S = struct('w', ones(nValid,1), 'robust_s', 1, 'covb', eye(nCond), 'resid', zeros(nValid,1));
        found_order = 0;  % AR order found by BIC on first iteration

        while norm(B - B0) / norm(B0) > 1e-2 && iter < maxiter
            B0 = B;
            res = y - X * B;

            % AR model fitting
            if maxBICSearch > 0 && ~nosearch
                if iter == 0
                    % First iteration: BIC search capped at maxBICSearch
                    a = nirs.math.ar_fit(res, min(p_i, maxBICSearch), false);
                    found_order = length(a) - 1;
                else
                    % Later iterations: direct fit at found order (skip BIC)
                    a = nirs.math.ar_fit(res, found_order, true);
                end
            else
                a = nirs.math.ar_fit(res, p_i, nosearch);
            end

            f = [1; -a(2:end)];

            % Filter design matrix and data
            Xf = myFilterLocal(f, X);

            if any(lstInvalid)
                yy = y;
                yy(lstInvalid) = interp1(find(lstValid), y(lstValid), find(lstInvalid), 'spline', true);
                yf = myFilterLocal(f, yy);
            else
                yf = myFilterLocal(f, y);
            end

            % Robust regression
            [B, S] = nirs.math.robustfit(Xf(lstValid, :), yf(lstValid), 'bisquare', tune, 'off');
            iter = iter + 1;
        end

        Xf_valid = Xf(lstValid, :);

        % Weighted design matrix
        w_col = S.w(:);  % ensure column
        wXf = bsxfun(@times, Xf_valid, w_col);

        % ============================================================
        % FAST Satterthwaite DOF (no n×n matrix)
        % ============================================================
        dfe_i = nirs.math.satterthwaite_dfe(w_col, wXf);

        % Store channel results
        dfe_all(i) = dfe_i;
        beta_all(:, i) = B;
        P_all(i) = length(a) - 1;

        L = pinv(Xf_valid' * Xf_valid);
        covb_per(:, :, i) = L * S.robust_s^2;

        Xfall_i = nan(nTime, size(wXf, 2));
        Xfall_i(lstValid, :) = wXf;
        Xfall{i} = Xfall_i;

        w_i = nan(nTime, 1);
        w_i(lstValid) = w_col;
        w_all(:, i) = w_i;
        a_all{i} = a;
        sigma2_all(i) = S.robust_s^2;

        tstat_all(:, i) = B ./ sqrt(diag(S.covb));
        pval_all(:, i) = 2 * tcdf(-abs(tstat_all(:, i)), dfe_i);
        ppos_all(:, i) = tcdf(-tstat_all(:, i), dfe_i);
        pneg_all(:, i) = tcdf(tstat_all(:, i), dfe_i);

        resid_i = nan(nTime, 1);
        resid_i(lstValid) = S.resid(:) .* w_col;
        resid(:, i) = resid_i;

        filter_all{i} = f;
        yf_valid = yf(lstValid);
        sse = norm(yf_valid - Xf_valid * B)^2;
        sst = norm(yf_valid - mean(yf_valid))^2;
        R2_all(i) = 1 - sse / sst;

        N = nValid - size(B, 2) - P_all(i);
        logLik_all(i) = (-N/2) * log(2*pi*sse/N) - N/2;
    end

    % ============================================================
    % CROSS-CHANNEL COVARIANCE (outside parfor)
    % ============================================================
    resid = bsxfun(@minus, resid, nanmedian(resid, 1));

    if nChan > 400
        C = resid' * resid;
    else
        C = zeros(nChan);
        for i = 1:nChan
            a_col = resid(:, i);
            for j = i:nChan
                C(i, j) = 1.4810 * nanmedian(a_col .* resid(:, j));
                C(j, i) = C(i, j);
            end
        end
    end
    C = C * (nanmean(sigma2_all ./ diag(C)));

    % Build 4D covariance tensor (batch Gram matrix, upper-triangle, backslash)
    covb = nirs.math.build_covb_tensor(Xfall, C);

    % ============================================================
    % Assemble output struct (compatible with original ar_irls)
    % ============================================================
    stats.beta = beta_all;
    stats.tstat = tstat_all;
    stats.pval = pval_all;
    stats.ppos = ppos_all;
    stats.pneg = pneg_all;
    stats.P = P_all;
    stats.w = w_all;
    stats.dfe = dfe_all;
    stats.covb = real(covb);
    stats.a = a_all;
    stats.sigma2 = sigma2_all;
    stats.filter = filter_all;
    stats.R2 = R2_all;
    stats.logLik = logLik_all;
end

function out = myFilterLocal(f, y)
    y1 = y(1, :);
    y = bsxfun(@minus, y, y1);
    out = filter(f, 1, y);
    out = bsxfun(@plus, out, sum(f) * y1);
end
