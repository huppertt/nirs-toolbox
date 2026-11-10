function test_bic_cap()
% TEST_BIC_CAP - Compare ar_irls_fast accuracy with different BIC cap settings
%
% Tests maxBICSearch = 0 (full), 3, 5, 8 and reports:
%   - Beta correlation with full search
%   - T-stat correlation with full search
%   - AR orders selected
%   - Timing

    fprintf('=== BIC Cap Accuracy Test ===\n\n');

    % Use simData for realistic fNIRS data
    raw = nirs.testing.simData;
    Fs = raw.Fs;
    Pmax = ceil(4 * Fs);
    fprintf('Fs = %.1f Hz, Pmax = %d\n\n', Fs, Pmax);

    d = raw.data;
    d = d - ones(size(d,1),1) * nanmean(d,1);

    % Build design matrix using public API (not protected methods)
    j = nirs.modules.AR_IRLS();
    [X, names] = nirs.design.createDesignMatrix(raw.stimulus, raw.time, j.basis);
    C = j.trend_func(raw.time);
    Xfull = [X C];

    caps = [0, 3, 5, 8];
    results = struct();

    for ci = 1:length(caps)
        cap = caps(ci);
        label = sprintf('maxBIC=%d', cap);
        if cap == 0, label = 'full (maxBIC=0)'; end

        t0 = tic;
        [stats, ~] = nirs.math.ar_irls_fast(d, Xfull, Pmax, [], false, false, false, cap);
        elapsed = toc(t0);

        results(ci).cap = cap;
        results(ci).label = label;
        results(ci).beta = stats.beta;
        results(ci).tstat = stats.tstat;
        results(ci).P = stats.P;
        results(ci).elapsed = elapsed;
        results(ci).dfe = stats.dfe;

        fprintf('%s:\n', label);
        fprintf('  Time:     %.3f s\n', elapsed);
        fprintf('  AR orders: min=%d  max=%d  mean=%.1f  median=%d\n', ...
            min(stats.P), max(stats.P), mean(stats.P), median(stats.P));
        fprintf('  Mean dfe:  %.1f\n', mean(stats.dfe));
        fprintf('\n');
    end

    % Compare all against full search (cap=0)
    fprintf('--- Accuracy vs full BIC search ---\n');
    ref_beta = results(1).beta(:);
    ref_tstat = results(1).tstat(:);
    ref_valid_b = ~isnan(ref_beta);
    ref_valid_t = ~isnan(ref_tstat);

    for ci = 2:length(caps)
        b = results(ci).beta(:);
        t = results(ci).tstat(:);

        beta_corr = corr(ref_beta(ref_valid_b), b(ref_valid_b));
        tstat_corr = corr(ref_tstat(ref_valid_t), t(ref_valid_t));
        beta_maxdiff = max(abs(ref_beta(ref_valid_b) - b(ref_valid_b)));
        tstat_maxdiff = max(abs(ref_tstat(ref_valid_t) - t(ref_valid_t)));
        speedup = results(1).elapsed / results(ci).elapsed;

        fprintf('%s vs full:\n', results(ci).label);
        fprintf('  Beta  corr: %.8f   max|diff|: %.6f\n', beta_corr, beta_maxdiff);
        fprintf('  Tstat corr: %.8f   max|diff|: %.6f\n', tstat_corr, tstat_maxdiff);
        fprintf('  Speedup:    %.2fx\n', speedup);
        fprintf('\n');
    end

    % Also compare against the ORIGINAL ar_irls
    fprintf('--- Accuracy vs original ar_irls ---\n');
    t0 = tic;
    [stats_orig, ~] = nirs.math.ar_irls(d, Xfull, Pmax);
    t_orig = toc(t0);
    fprintf('Original ar_irls: %.3f s\n\n', t_orig);

    orig_beta = stats_orig.beta(:);
    orig_tstat = stats_orig.tstat(:);
    orig_valid_b = ~isnan(orig_beta);
    orig_valid_t = ~isnan(orig_tstat);

    for ci = 1:length(caps)
        b = results(ci).beta(:);
        t = results(ci).tstat(:);

        beta_corr = corr(orig_beta(orig_valid_b), b(orig_valid_b));
        tstat_corr = corr(orig_tstat(orig_valid_t), t(orig_valid_t));
        speedup = t_orig / results(ci).elapsed;

        fprintf('%s vs original:\n', results(ci).label);
        fprintf('  Beta  corr: %.8f\n', beta_corr);
        fprintf('  Tstat corr: %.8f\n', tstat_corr);
        fprintf('  Speedup:    %.2fx\n', speedup);
        fprintf('\n');
    end

end
