function stats = ar_irls_REML_fast(d, X, Pmax, tune, useGPU, singlePrecision, maxBICSearch)
% AR_IRLS_REML_FAST  Parallelized AR-IRLS with REML robust fitting and fast
%   Satterthwaite degrees-of-freedom estimation.
%
%   stats = ar_irls_REML_fast(d, X, Pmax, tune, useGPU, singlePrecision, maxBICSearch)
%
%   This is a drop-in replacement for ar_irls_REML that uses parfor for
%   channel-level parallelism and nirs.math.satterthwaite_dfe for O(n*k)
%   DOF computation instead of forming the full n-by-n hat matrix.
%
%   Inputs:
%       d               - nTime x nChan data matrix
%       X               - nTime x nCond design matrix
%       Pmax            - maximum AR model order (scalar or per-channel vector)
%       tune            - tuning constant for bisquare robust weighting (default 4.685)
%       useGPU          - unused, kept for signature compatibility (default false)
%       singlePrecision - unused, kept for signature compatibility (default false)
%       maxBICSearch    - max AR orders for BIC search (0=full, default 0)
%
%   Output:
%       stats - struct with fields: beta, tstat, pval, ppos, pneg, P, w,
%               dfe, covb, a, sigma2, filter, R2

    warning('off', 'stats:statrobustfit:IterationLimit');
    if nargin < 4 || isempty(tune), tune = 4.685; end
    if nargin < 5, useGPU = false; end   %#ok<NASGU> kept for API compat
    if nargin < 6, singlePrecision = false; end   %#ok<NASGU>
    if nargin < 7 || isempty(maxBICSearch), maxBICSearch = 0; end

    nCond = size(X, 2);
    nChan = size(d, 2);
    nTime = size(d, 1);

    % ------------------------------------------------------------------
    %  Pre-allocate temporary arrays (parfor-friendly, no struct indexing)
    % ------------------------------------------------------------------
    beta_all    = zeros(nCond, nChan);
    tstat_all   = zeros(nCond, nChan);
    pval_all    = zeros(nCond, nChan);
    ppos_all    = zeros(nCond, nChan);
    pneg_all    = zeros(nCond, nChan);
    P_all       = zeros(nChan, 1);
    w_all       = zeros(nTime, nChan);
    dfe_all     = zeros(nChan, 1);
    covb_per    = zeros(nCond, nCond, nChan);
    sigma2_all  = zeros(nChan, 1);
    a_all       = cell(nChan, 1);
    filter_all  = cell(nChan, 1);
    R2_all      = zeros(nChan, 1);
    resid       = nan(nTime, nChan);
    Xfall       = cell(nChan, 1);

    % Pmax: allow scalar or per-channel vector
    if length(Pmax) == 1
        Pmax_vec = repmat(Pmax, nChan, 1);
    else
        Pmax_vec = Pmax(:);
    end

    % ------------------------------------------------------------------
    %  Channel loop (parfor)
    % ------------------------------------------------------------------
    parfor i = 1:nChan
        y = d(:, i);

        % --- skip all-zero channels ---
        if ~any(y)
            Xfall{i}      = nan(nTime, nCond);
            beta_all(:, i) = NaN;
            covb_per(:, :, i) = NaN;
            w_all(:, i)   = NaN;
            a_all{i}      = [];
            sigma2_all(i) = NaN;
            tstat_all(:, i) = NaN;
            pval_all(:, i)  = NaN;
            ppos_all(:, i)  = NaN;
            pneg_all(:, i)  = NaN;
            resid(:, i)     = NaN;
            filter_all{i}   = [];
            R2_all(i)       = NaN;
            dfe_all(i)      = NaN;
            continue;
        end

        lstValid   = ~isnan(y);
        lstInvalid = ~lstValid;

        % Initial OLS fit (backslash, matching original REML)
        B  = X(lstValid, :) \ y(lstValid);
        B0 = 1e6 * ones(size(B));

        p_i     = Pmax_vec(i);
        iter    = 0;
        maxiter = 10;
        f  = 1;
        Xf = X;
        yf = y;
        S  = struct('w', ones(sum(lstValid), 1), 'covb', eye(nCond), ...
                    'sigma', 1, 'resid', zeros(sum(lstValid), 1));
        found_order = 0;

        while norm(B - B0) / norm(B0) > 1e-2 && iter < maxiter
            B0  = B;
            res = y - X * B;

            % AR fit with optional BIC cap + skip-after-first
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
            f = [1; -a(2:end)];

            Xf = myFilterLocal(f, X);

            if any(lstInvalid)
                yy = y;
                yy(lstInvalid) = interp1(find(lstValid), y(lstValid), ...
                                          find(lstInvalid), 'spline', true);
                yf = myFilterLocal(f, yy);
            else
                yf = myFilterLocal(f, y);
            end

            % REML robust fit
            [B, S] = nirs.math.robustfit_reml(Xf(lstValid, :), yf(lstValid), ...
                                              'bisquare', tune, 'off');
            iter = iter + 1;
        end

        % --- Post-convergence quantities ---
        Xf_valid = Xf(lstValid, :);
        w_col    = S.w(:);
        wXf      = bsxfun(@times, Xf_valid, w_col);

        % Fast Satterthwaite DOF (O(n*k) instead of O(n^2))
        dfe_i      = nirs.math.satterthwaite_dfe(w_col, wXf);
        dfe_all(i) = dfe_i;

        % Store per-channel results
        beta_all(:, i) = B;
        P_all(i)       = length(a) - 1;

        % covb: store S.covb directly (REML-specific)
        covb_per(:, :, i) = S.covb;

        Xfall_i = nan(nTime, size(wXf, 2));
        Xfall_i(lstValid, :) = wXf;
        Xfall{i} = Xfall_i;

        w_i = nan(nTime, 1);
        w_i(lstValid) = w_col;
        w_all(:, i) = w_i;

        a_all{i}      = a;
        sigma2_all(i) = S.sigma^2;   % REML uses S.sigma

        tstat_all(:, i) = B ./ sqrt(diag(S.covb));
        pval_all(:, i)  = 2 * tcdf(-abs(tstat_all(:, i)), dfe_i);
        ppos_all(:, i)  = tcdf(-tstat_all(:, i), dfe_i);
        pneg_all(:, i)  = tcdf(tstat_all(:, i), dfe_i);

        resid_i = nan(nTime, 1);
        resid_i(lstValid) = S.resid(:) .* w_col;
        resid(:, i) = resid_i;

        filter_all{i} = f;

        % R2: MAD-based formula (REML-specific)
        yf_valid = yf(lstValid);
        R2_all(i) = max(1 - mad(yf_valid - Xf_valid * B) / mad(yf_valid), 0);
    end

    % ------------------------------------------------------------------
    %  Aggregate scalar dfe (mean over valid channels)
    % ------------------------------------------------------------------
    dfe_scalar = nanmean(dfe_all);

    % ------------------------------------------------------------------
    %  Pack per-channel stats into output struct
    % ------------------------------------------------------------------
    stats.beta    = beta_all;
    stats.tstat   = tstat_all;
    stats.pval    = pval_all;
    stats.ppos    = ppos_all;
    stats.pneg    = pneg_all;
    stats.P       = P_all;
    stats.w       = w_all;
    stats.dfe     = dfe_scalar;
    stats.covb    = covb_per;       % will be overwritten below with 4-D version
    stats.a       = a_all;
    stats.sigma2  = sigma2_all(:)';
    stats.filter  = filter_all;
    stats.R2      = R2_all(:)';

    % ------------------------------------------------------------------
    %  Cross-channel covariance using ar_corr
    % ------------------------------------------------------------------
    covb4 = zeros(nCond, nCond, nChan, nChan);
    C = nirs.sFC.ar_corr(resid, Pmax, true);

    for i = 1:nChan
        for j = i:nChan
            covb4(:, :, i, j) = covb4(:, :, i, j) + ...
                C(i, j) * sqrt(covb_per(:, :, i) .* covb_per(:, :, j));
        end
    end
    stats.covb = real(covb4);

    % ------------------------------------------------------------------
    %  Re-compute tstat / pval from 4-D covb (upper-triangle fill)
    % ------------------------------------------------------------------
    for i = 1:nChan
        stats.tstat(:, i) = stats.beta(:, i) ./ ...
            sqrt(diag(squeeze(stats.covb(:, :, i, i))));
        stats.pval(:, i)  = 2 * tcdf(-abs(stats.tstat(:, i)), stats.dfe);
        stats.ppos(:, i)  = tcdf(-stats.tstat(:, i), stats.dfe);
        stats.pneg(:, i)  = tcdf(stats.tstat(:, i), stats.dfe);
    end
end

% ======================================================================
%  Local filter helper (must be defined inside file for parfor visibility)
% ======================================================================
function out = myFilterLocal(f, y)
    y1  = y(1, :);
    y   = bsxfun(@minus, y, y1);
    out = filter(f, 1, y);
    out = bsxfun(@plus, out, sum(f) * y1);
end
