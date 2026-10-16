function stats = ar_irls_Mixed( d,X,Z,Pmax,tune,maxBICSearch)
% See the following for the related publication:
% http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3756568/
%
% d is matrix containing the data; each column is a channel of data
%
% X is the regression/design matrix
%
% Pmax is the maximum AR model order that you want to consider. A
% purely speculative guess is that the average model order is
% approximatley equal to the 2-3 times the sampling rate, so setting Pmax
% to 4 or 5 times the sampling rate should work fine.  The code does not
% suffer a hugeperformance hit by using a higher Pmax; however, the number
% of time points used to estimate the AR model will be
% "# of time points - Pmax", so don't set Pmax too high.
%
% "tune" is the tuning constant used for Tukey's bisquare function during
% iterative reweighted least squares. The default value is 4.685.
% Decreasing "tune" will make the regression less sensitive to outliers,
% but at the expense of performance (statistical efficiency) when data
% does not have outliers. For reference, the values of tune for 85, 90,
% and 95% statistical efficiency are
%
% tune = 4.685 --> 95%
% tune = 4.00  --> ~90%
% tune = 3.55  --> ~85%
%
% I have not tested these to find an optimal value for the "average" NIRS
% dataset; however, 4.685 was used in the published simulations and worked
% quite well even with a high degree of motion artifacts from children.
% If you really want to adjust it, you could use the above values as a
% guideline.
%
% DO NOT preprocess your data with a low pass filter.
% The algorithm is trying to transform the residual to create a
% white spectrum.  If part of the spectrum is missing due to low pass
% filtering, the AR coefficients will be unstable.  High pass filtering
% may be ok, but I suggest putting orthogonal polynomials (e.g. Legendre)
% or low frequency discrete cosine terms directly into the design matrix
% (e.g. from spm_dctmtx.m from SPM).  Don't use regular polynomials
% (e.g. 1 t t^2 t^3 t^4 ...) as this can result in a poorly conditioned
% design matrix.
%
% If you choose to resample your data to a lower sampling frequency,
% makes sure to choose an appropriate cutoff frequency so that that the
% resulting time series is not missing part of the frequency spectrum
% (up to the Nyquist bandwidth).  The code should work fine on 10-30 Hz
% data.

warning('off','stats:statrobustfit:IterationLimit')

if nargin < 5
    tune = 4.685;
end
if nargin < 6 || isempty(maxBICSearch), maxBICSearch = 0;
end

% preallocate stats
nCond = size(X,2);
nChan = size(d,2);
nTime = size(d,1);

[U,S,V]=nirs.math.mysvd(X);
lstgood = find(diag(S)>eps(1)*50);
X=U(:,lstgood)*S(lstgood,lstgood);

nCondOut = size(V,2);  % full original condition count after SVD projection

% pre-allocate sliced output arrays for parfor compatibility
beta_all   = zeros(nCondOut, nChan);
tstat_all  = zeros(nCondOut, nChan);
pval_all   = zeros(nCondOut, nChan);
ppos_all   = zeros(nCondOut, nChan);
pneg_all   = zeros(nCondOut, nChan);
covb_all   = zeros(nCondOut, nCondOut, nChan);
dfe_all    = zeros(1, nChan);
P_all      = zeros(1, nChan);
sigma2_all = zeros(1, nChan);
R2_all     = zeros(1, nChan);
w_all      = zeros(nTime, nChan);
a_all      = cell(1, nChan);
filter_all = cell(1, nChan);

% loop through each channel (parallelized)
parfor i = 1:nChan
    y = d(:,i);

    % initial fit
    B = X \ y;
    B0 = 1e6*ones(size(B));
    LME = struct('residuals', y - X*B);

    % iterative re-weighted least squares
    iter = 0;
    maxiter = 5;

    a = [];
    f = [];
    w = [];
    yf = [];
    found_order = 0;

    % while our coefficients are changing greater than some threshold
    % and it's less than the max number of iterations
    while norm(B-B0)/norm(B0) > 1e-2 && iter < maxiter
        % store the last fit
        B0 = B;

        % get the residual
        res = LME.residuals;

        % fit the residual to an ar model
        if maxBICSearch > 0
            if iter == 0
                a = nirs.math.ar_fit(res, min(Pmax, maxBICSearch), false);
                found_order = length(a) - 1;
            else
                a = nirs.math.ar_fit(res, found_order, true);
            end
        else
            a = nirs.math.ar_fit(res, Pmax);
        end

        % create a whitening filter from the coefficients
        f = [1; -a(2:end)];

        % filter the design matrix
        Xf = myFilter(f,X);
        if(~isempty(Z))
            Zf = myFilter(f,Z);
        else
            Zf=[];
        end

        % subtract constant from AR model and filter the data
        yf = myFilter(f,y);

        for id=1:4
            w = wfun(LME.residuals,tune);
            LME = fitlmematrix(Xf,yf,Zf,[],'FitMethod','REML','weights',w);
        end
        B=LME.Coefficients.Estimate;
        iter = iter + 1;
    end

    %  Satterthwaite estimate of model DOF
    dfe_i = min(LME.Coefficients.DF);

    % moco data & statistics
    beta_i = V(lstgood,:)'*LME.Coefficients.Estimate;
    covb_i = V(lstgood,:)'*LME.CoefficientCovariance*V(lstgood,:);

    tstat_i = beta_i ./ diag(sqrt(covb_i));
    pval_i  = 2*tcdf(-abs(tstat_i), dfe_i);        % two-sided
    ppos_i  = tcdf(-tstat_i, dfe_i);                % one-sided (positive only)
    pneg_i  = tcdf(tstat_i, dfe_i);                 % one-sided (negative only)

    % store into sliced arrays
    beta_all(:,i)    = beta_i;
    tstat_all(:,i)   = tstat_i;
    pval_all(:,i)    = pval_i;
    ppos_all(:,i)    = ppos_i;
    pneg_all(:,i)    = pneg_i;
    covb_all(:,:,i)  = covb_i;
    dfe_all(i)       = dfe_i;
    P_all(i)         = length(a)-1;
    sigma2_all(i)    = mad(LME.residuals)^2;
    R2_all(i)        = max(1-mad(yf-LME.residuals)/mad(yf),0);
    w_all(:,i)       = w;
    a_all{i}         = a;
    filter_all{i}    = f;
end

% assemble stats struct from parallel results
stats.beta   = beta_all;
stats.tstat  = tstat_all;
stats.pval   = pval_all;
stats.ppos   = ppos_all;
stats.pneg   = pneg_all;
stats.covb   = covb_all;
stats.dfe    = dfe_all;
stats.P      = P_all;
stats.sigma2 = sigma2_all;
stats.R2     = R2_all;
stats.w      = w_all;
stats.a      = a_all;
stats.filter = filter_all;
end

%%
function out = myFilter( f, y )
% here we are just making the first value zero before filtering to
% avoid weird effects introduced by zero padding
y1 = y(1,:);

y = bsxfun(@minus,y,y1);

out = filter(f, 1, y);
out = bsxfun(@plus,out,sum(f)*y1); % add the corrected offset back

end



function w = wfun(r,tune)
s = mad(r, 0) / 0.6745;
r = r/s/tune;

w = (1 - r.^2) .* (r < 1 & r > -1);
end