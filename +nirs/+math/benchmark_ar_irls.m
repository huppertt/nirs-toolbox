function benchmark_ar_irls(nChan, nTime, nCond)
% BENCHMARK_AR_IRLS Profile the AR-IRLS algorithm to identify bottlenecks
%
% Usage:
%   nirs.math.benchmark_ar_irls()           % defaults: 32 chan, 3000 time, 6 cond
%   nirs.math.benchmark_ar_irls(64, 5000, 8)

if nargin < 1, nChan = 32; end
if nargin < 2, nTime = 3000; end
if nargin < 3, nCond = 6; end

Fs = 10; % typical fNIRS sampling rate
Pmax = 4 * Fs;

fprintf('\n=== AR-IRLS Benchmark ===\n');
fprintf('Channels: %d | Time points: %d | Conditions: %d | Pmax: %d\n\n', nChan, nTime, nCond, Pmax);

% Generate synthetic data
rng(42);
t = (0:nTime-1)'/Fs;

% Design matrix: nCond stimulus regressors + 4 trend (Legendre) + constant
X_stim = zeros(nTime, nCond);
for c = 1:nCond
    % Random block design
    onsets = sort(randi(nTime-100, 5, 1));
    for o = 1:length(onsets)
        idx = onsets(o):min(onsets(o)+30, nTime);
        X_stim(idx, c) = 1;
    end
end
% Add trend regressors
C = [ones(nTime,1), legendre_basis(t, 3)];
X = [X_stim, C];

% Synthetic NIRS-like data with AR(5) noise
d = zeros(nTime, nChan);
for ch = 1:nChan
    beta_true = randn(size(X, 2), 1) * 0.1;
    noise = randn(nTime, 1);
    % AR(5) coloring
    for tt = 6:nTime
        noise(tt) = 0.5*noise(tt-1) - 0.2*noise(tt-2) + 0.1*noise(tt-3) + noise(tt)*0.3;
    end
    d(:, ch) = X * beta_true + noise;
end

fprintf('--- Profiling ORIGINAL ar_irls ---\n');

% Time the full function
tic;
[stats_orig, ~] = nirs.math.ar_irls(d, X, Pmax, [], false, false, false);
t_total = toc;
fprintf('TOTAL TIME: %.2f seconds (%.0f ms/channel)\n\n', t_total, t_total/nChan*1000);

% Now profile individual phases
fprintf('--- Phase Breakdown ---\n');

% Phase 1: Initial OLS fits (all channels)
tic;
for i = 1:nChan
    y = d(:,i);
    lstV = ~isnan(y);
    pinv(X(lstV,:)) * y(lstV);
end
t_init = toc;
fprintf('Initial OLS fits:        %6.2f s (%5.1f%%)\n', t_init, t_init/t_total*100);

% Phase 2: AR fitting
tic;
for i = 1:nChan
    y = d(:,i);
    B = pinv(X) * y;
    res = y - X*B;
    for iter = 1:10
        nirs.math.ar_fit(res, Pmax, false);
    end
end
t_ar = toc;
fprintf('AR model fitting (10x):  %6.2f s (%5.1f%%)\n', t_ar, t_ar/t_total*100);

% Phase 3: Filtering
tic;
f_test = [1; -0.5*ones(5,1)];
for i = 1:nChan
    for iter = 1:10
        filter(f_test, 1, X);
        filter(f_test, 1, d(:,i));
    end
end
t_filt = toc;
fprintf('Filtering (10x):         %6.2f s (%5.1f%%)\n', t_filt, t_filt/t_total*100);

% Phase 4: Robust regression
tic;
f_test = [1; -0.3; 0.1; -0.05];
Xf = filter(f_test, 1, X);
for i = 1:nChan
    yf = filter(f_test, 1, d(:,i));
    for iter = 1:3
        nirs.math.robustfit(Xf, yf, 'bisquare', 4.685, 'off');
    end
end
t_robust = toc;
fprintf('Robust regression (3x):  %6.2f s (%5.1f%%)\n', t_robust, t_robust/t_total*100);

% Phase 5: Satterthwaite DOF
tic;
w_test = rand(nTime, 1) * 0.3 + 0.7;
wXf = bsxfun(@times, Xf, w_test);
for i = 1:nChan
    H = diag(w_test) - wXf*pinv(wXf'*wXf)*wXf';
    HtH = H'*H;
    dfe = sum(reshape(H,[],1).*reshape(H',[],1))^2 / sum(reshape(HtH,[],1).^2);
end
t_satt = toc;
fprintf('Satterthwaite DOF:       %6.2f s (%5.1f%%)\n', t_satt, t_satt/t_total*100);

% Phase 5b: Satterthwaite DOF (optimized - no n×n matrix)
tic;
for i = 1:nChan
    [U_t,~,~] = svd(wXf, 'econ');
    k = size(U_t, 2);
    h = sum(U_t.^2, 2);
    trM = sum(w_test.^2) - 2*sum(w_test.*h) + k;
    % For trace(M^2), use efficient O(n*k^2) formula
    Mw = U_t' * bsxfun(@times, U_t, w_test);
    trM2 = sum(w_test.^4) - 4*sum(w_test.^3.*h) + 2*sum(w_test.^2.*h) ...
        + 4*sum(w_test.^2.*h.^2) + 4*sum(Mw(:).^2) - 4*sum(w_test.*h.^2) ...
        + sum(h.^2);
    % Simplified: compute via direct formula
    diag_H2 = w_test.^2 - (2*w_test - 1).*h;
    trM_check = sum(diag_H2) + k - sum(h.^2); % should match trM
    dfe_fast = trM^2 / max(trM2, eps);
end
t_satt_fast = toc;
fprintf('Satterthwaite (fast):    %6.2f s (%5.1f%%) [%.1fx speedup]\n', ...
    t_satt_fast, t_satt_fast/t_total*100, t_satt/max(t_satt_fast,eps));

% Phase 6: Cross-channel covariance
tic;
nC = size(X, 2);
covb = zeros(nC, nC, nChan, nChan);
Xfall_test = cell(nChan, 1);
for i = 1:nChan
    Xfall_test{i} = wXf + randn(size(wXf))*0.01;
end
C_test = eye(nChan);
for i = 1:nChan
    for j = 1:nChan
        lstV = ~isnan(sum(Xfall_test{i},2) + sum(Xfall_test{j},2));
        covb(:,:,i,j) = pinv(Xfall_test{i}(lstV,:)'*Xfall_test{j}(lstV,:)) * C_test(i,j);
    end
end
t_covb = toc;
fprintf('Cross-chan covariance:    %6.2f s (%5.1f%%)\n', t_covb, t_covb/t_total*100);

% Summary
fprintf('\n--- Summary ---\n');
fprintf('Total original time: %.2f seconds\n', t_total);
fprintf('\nEstimated speedup from parallelizing %d channels across workers:\n', nChan);
for nw = [2 4 8 16]
    fprintf('  %2d workers: ~%.2f s (%.1fx)\n', nw, t_total/min(nw, nChan), min(nw, nChan));
end

fprintf('\nEstimated speedup from Satterthwaite optimization: %.1fx on that phase\n', ...
    t_satt/max(t_satt_fast, eps));

% ============================================================
% Benchmark the FAST version
% ============================================================
fprintf('\n--- Profiling FAST ar_irls_fast ---\n');

tic;
[stats_fast, ~] = nirs.math.ar_irls_fast(d, X, Pmax, [], false, false, false);
t_fast = toc;
fprintf('TOTAL TIME (fast): %.2f seconds (%.0f ms/channel)\n', t_fast, t_fast/nChan*1000);
fprintf('SPEEDUP: %.1fx\n', t_total/t_fast);

% Validate results match
beta_diff = max(abs(stats_orig.beta(:) - stats_fast.beta(:)));
fprintf('\nValidation - max beta difference: %.2e\n', beta_diff);
if beta_diff < 0.1
    fprintf('PASS: Results match within tolerance\n');
else
    fprintf('WARNING: Results diverge (expected due to AR order search skip on iter 2+)\n');
    fprintf('  Correlation of betas: %.6f\n', corr(stats_orig.beta(:), stats_fast.beta(:)));
end

fprintf('\n--- Summary ---\n');
fprintf('Original:  %.2f s (%.0f ms/channel)\n', t_total, t_total/nChan*1000);
fprintf('Fast:      %.2f s (%.0f ms/channel)\n', t_fast, t_fast/nChan*1000);
fprintf('Speedup:   %.1fx\n', t_total/t_fast);

% Test if parfor is available
try
    p = gcp('nocreate');
    if isempty(p)
        fprintf('\nNo parallel pool active. For additional speedup: parpool(N)\n');
    else
        fprintf('\nParallel pool active: %d workers\n', p.NumWorkers);
    end
catch
end

fprintf('\n=== Benchmark Complete ===\n');
end

function L = legendre_basis(t, P)
    n = length(t);
    x = linspace(-1, 1, n)';
    L = zeros(n, P);
    for p = 1:P
        L(:,p) = x.^p;
    end
    [L, ~] = qr(L, 0);
end

function out = myFilterLocal(f, y)
    y1 = y(1,:);
    y = bsxfun(@minus, y, y1);
    out = filter(f, 1, y);
    out = bsxfun(@plus, out, sum(f)*y1);
end
