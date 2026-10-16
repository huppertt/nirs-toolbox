function stats = ar_irnsls_fast(d, X, Pmax, tune, useGPU, singlePrecision, maxBICSearch)
% AR_IRNSLS_FAST - Parallelized AR(P)-IRNSLS regression
%
% Drop-in replacement for nirs.math.ar_irnsls with these optimizations:
%   1. parfor across channels (falls back to for-loop if no pool)
%   2. Satterthwaite DOF without forming n*n H matrix (O(n*k^2) vs O(n^2))
%   3. Diagonal-only cross-channel covariance (matching original)
%   4. Capped BIC search + skip BIC after first IRLS iteration
%
% Usage:
%   stats = nirs.math.ar_irnsls_fast(d, X, Pmax, tune, useGPU, singlePrecision, maxBICSearch)
%
% See also: nirs.math.ar_irnsls, nirs.math.satterthwaite_dfe

    warning('off', 'stats:statrobustfit:IterationLimit');

    if nargin < 4 || isempty(tune), tune = 4.685; end
    if nargin < 5, useGPU = false; end
    if nargin < 6, singlePrecision = false; end
    if nargin < 7 || isempty(maxBICSearch), maxBICSearch = 0; end

    nCond = size(X, 2);
    nChan = size(d, 2);
    nTime = size(d, 1);

    % Preallocate outputs as simple arrays for parfor compatibility
    beta_all   = zeros(nCond, nChan);
    tstat_all  = zeros(nCond, nChan);
    pval_all   = zeros(nCond, nChan);
    ppos_all   = zeros(nCond, nChan);
    pneg_all   = zeros(nCond, nChan);
    P_all      = zeros(nChan, 1);
    w_all      = zeros(nTime, nChan);
    dfe_all    = zeros(nChan, 1);
    covb_per   = zeros(nCond, nCond, nChan);
    sigma2_all = zeros(nChan, 1);
    a_all      = cell(nChan, 1);
    filter_all = cell(nChan, 1);
    R2_all     = zeros(nChan, 1);

    % Handle per-channel Pmax (scalar only expected, but vec for safety)
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
            beta_all(:, i)    = NaN;
            covb_per(:, :, i) = NaN;
            w_all(:, i)       = NaN;
            a_all{i}          = [];
            sigma2_all(i)     = NaN;
            tstat_all(:, i)   = NaN;
            pval_all(:, i)    = NaN;
            ppos_all(:, i)    = NaN;
            pneg_all(:, i)    = NaN;
            filter_all{i}     = [];
            R2_all(i)         = NaN;
            continue;
        end

        lstValid   = ~isnan(y);
        lstInvalid = ~lstValid;
        p_i = Pmax_vec(i);

        % Initial OLS fit (backslash, matching original)
        B  = X(lstValid, :) \ y(lstValid);
        B0 = 1e6 * ones(size(B));

        iter    = 0;
        maxiter = 10;
        f  = 1;
        Xf = X;
        yf = y;

        % Default S struct for parfor (nonstationary_robustfit fields)
        nValid = sum(lstValid);
        S = struct('w', ones(nValid, 1), ...
                   'robust_s', 1, ...
                   'covb', eye(nCond), ...
                   'resid', zeros(nValid, 1), ...
                   'dfe', nValid - nCond);

        found_order = 0;

        while norm(B - B0) / norm(B0) > 1e-2 && iter < maxiter
            B0 = B;

            % Residual
            res = y - X * B;

            % AR model fitting with optional BIC cap + skip-after-first
            if maxBICSearch > 0
                if iter == 0
                    a = nirs.math.ar_fit(res, min(p_i, maxBICSearch), false);
                    found_order = length(a) - 1;
                else
                    a = nirs.math.ar_fit(res, found_order, true);
                end
            else
                a = nirs.math.ar_fit(res, p_i);
            end

            % Whitening filter from AR coefficients
            f = [1; -a(2:end)];

            % Filter design matrix
            Xf = myFilterLocal(f, X);

            % Filter data (interpolate over NaNs if needed)
            if any(lstInvalid)
                yy = y;
                yy(lstInvalid) = interp1(find(lstValid), y(lstValid), ...
                                         find(lstInvalid), 'spline', true);
                yf = myFilterLocal(f, yy);
            else
                yf = myFilterLocal(f, y);
            end

            % Nonstationary robust regression
            [B, S] = nirs.math.nonstationary_robustfit( ...
                         Xf(lstValid, :), yf(lstValid), ...
                         'bisquare', tune, 'off');

            iter = iter + 1;
        end

        Xf_valid = Xf(lstValid, :);
        yf_valid = yf(lstValid);

        % Weighted design matrix
        w_col = S.w(:);
        wXf = bsxfun(@times, Xf_valid, w_col);

        % ============================================================
        % FAST Satterthwaite DOF (no n*n matrix)
        % ============================================================
        dfe_i = nirs.math.satterthwaite_dfe(w_col, wXf);

        % Store channel results
        dfe_all(i)        = dfe_i;
        beta_all(:, i)    = B;
        P_all(i)          = length(a) - 1;
        covb_per(:, :, i) = S.covb;
        sigma2_all(i)     = S.robust_s^2;

        w_i = nan(nTime, 1);
        w_i(lstValid) = w_col;
        w_all(:, i) = w_i;

        a_all{i}      = a;
        filter_all{i} = f;

        tstat_all(:, i) = B ./ sqrt(diag(S.covb));
        pval_all(:, i)  = 2 * tcdf(-abs(tstat_all(:, i)), dfe_i);
        ppos_all(:, i)  = tcdf(-tstat_all(:, i), dfe_i);
        pneg_all(:, i)  = tcdf(tstat_all(:, i), dfe_i);

        R2_all(i) = max(1 - mad(yf_valid - Xf_valid * B) / mad(yf_valid), 0);
    end

    % ============================================================
    % CROSS-CHANNEL COVARIANCE - DIAGONAL ONLY (matching original)
    % ============================================================
    covb = zeros(nCond, nCond, nChan, nChan);
    for i = 1:nChan
        covb(:, :, i, i) = covb_per(:, :, i);
    end

    % ============================================================
    % Assemble output struct (compatible with original ar_irnsls)
    % ============================================================
    stats.beta   = beta_all;
    stats.tstat  = tstat_all;
    stats.pval   = pval_all;
    stats.ppos   = ppos_all;
    stats.pneg   = pneg_all;
    stats.P      = P_all;
    stats.w      = w_all;
    stats.dfe    = nanmean(dfe_all);
    stats.covb   = real(covb);
    stats.a      = a_all;
    stats.sigma2 = sigma2_all;
    stats.filter = filter_all;
    stats.R2     = R2_all;
end

function out = myFilterLocal(f, y)
    y1  = y(1, :);
    y   = bsxfun(@minus, y, y1);
    out = filter(f, 1, y);
    out = bsxfun(@plus, out, sum(f) * y1);
end
