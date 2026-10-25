function stats = ar_irls_priors( d,X,isRE,Pmax,tune,type,maxBICSearch )
% isFE=[1 0]; -  boolean array of if Fixed Effects (or random effects) in
% the model- defines which regressors have spatial covariance priors
% applied to them

warning('off','stats:statrobustfit:IterationLimit')

if nargin < 5
    tune = 4.685;
end
if nargin < 7 || isempty(maxBICSearch), maxBICSearch = 0; end

Priors=zeros(size(X,2),1);
Lambda=zeros(size(X,2));

nChan=size(d,2);
nCond=size(X,2);
nTime=size(d,1);

% Pre-allocate output arrays
beta_all   = zeros(nCond, nChan);
tstat_all  = zeros(nCond, nChan);
pval_all   = zeros(nCond, nChan);
ppos_all   = zeros(nCond, nChan);
pneg_all   = zeros(nCond, nChan);
P_all      = zeros(1, nChan);
covb_all   = zeros(nCond, nCond, nChan);
w_all      = zeros(nTime, nChan);
a_all      = cell(1, nChan);
sigma2_all = zeros(1, nChan);
filter_all = cell(1, nChan);
R2_all     = zeros(1, nChan);
dfe_all    = zeros(1, nChan);
Ball       = zeros(nCond, nChan);

utype=unique(type);
ntypes=length(utype);
disp('*');

for it=1:ntypes
    lst=find(ismember(type,utype{it}));

    Lambda0=1E4*eye(size(X,2));
    outeriter=0;

    while(norm(Lambda0-Lambda)>1E-4 & outeriter<5)
        Lambda0=Lambda;

        % Capture Priors and Lambda as local variables for parfor broadcast
        Priors_local = Priors;
        Lambda_local = Lambda;

        nLst = length(lst);

        % Temporary sliced arrays for parfor
        beta_tmp   = zeros(nCond, nLst);
        tstat_tmp  = zeros(nCond, nLst);
        pval_tmp   = zeros(nCond, nLst);
        ppos_tmp   = zeros(nCond, nLst);
        pneg_tmp   = zeros(nCond, nLst);
        P_tmp      = zeros(1, nLst);
        covb_tmp   = zeros(nCond, nCond, nLst);
        w_tmp      = zeros(nTime, nLst);
        a_tmp      = cell(1, nLst);
        sigma2_tmp = zeros(1, nLst);
        filter_tmp = cell(1, nLst);
        R2_tmp     = zeros(1, nLst);
        dfe_tmp    = zeros(1, nLst);
        Ball_tmp   = zeros(nCond, nLst);

        parfor ii = 1:nLst
            i = lst(ii);

            fprintf( 'Finished Iter-%4i %4i of %4i.\n', outeriter+1, i, nChan );

            y = d(:,i) - X*Priors_local;
            B = pinv(X'*X + Lambda_local) * X'*y;

            w_loc = zeros(nTime, 1);
            covb_loc = zeros(nCond, nCond);
            a_loc = [];
            f_loc = [];

            found_order = 0;
            for iter2=1:5
                res = y - X*B;
                if maxBICSearch > 0
                    if iter2 == 1
                        a_loc = nirs.math.ar_fit(res, min(Pmax, maxBICSearch), false);
                        found_order = length(a_loc) - 1;
                    else
                        a_loc = nirs.math.ar_fit(res, found_order, true);
                    end
                else
                    a_loc = nirs.math.ar_fit(res, Pmax);
                end
                f_loc = [1; -a_loc(2:end)];

                % filter the design matrix
                Xf = myFilterLocal(f_loc, X);

                % subtract constant from AR model and filter the data
                yf = myFilterLocal(f_loc, y);

                for iter=1:10
                    r = yf - Xf*B;
                    s = mad(r, 0) / 0.6745;
                    r = r/s/tune;
                    w_loc = (1 - r.^2) .* (r < 1 & r > -1);
                    Xfw = bsxfun(@times, Xf, w_loc);
                    yfw = bsxfun(@times, yf, w_loc);
                    B = pinv(Xfw'*Xfw + Lambda_local) * Xfw'*yfw;
                end

                r = yf - Xf*B;
                covb_loc = pinv(Xfw'*Xfw + Lambda_local) * mean(r.^2);
            end

            % Compute per-channel statistics
            tstat_loc = B ./ sqrt(diag(covb_loc));
            tstat_loc = real(tstat_loc);
            tstat_loc(isnan(tstat_loc)) = 0;
            dfe_loc   = length(yf) - length(B);

            beta_tmp(:,ii)     = B + Priors_local;
            tstat_tmp(:,ii)    = tstat_loc;
            pval_tmp(:,ii)     = 2*tcdf(-abs(tstat_loc), dfe_loc);
            ppos_tmp(:,ii)     = tcdf(-tstat_loc, dfe_loc);
            pneg_tmp(:,ii)     = tcdf(tstat_loc, dfe_loc);
            P_tmp(ii)          = length(a_loc) - 1;
            covb_tmp(:,:,ii)   = covb_loc;
            w_tmp(:,ii)        = w_loc;
            a_tmp{ii}          = a_loc;
            sigma2_tmp(ii)     = mad(r)^2;
            filter_tmp{ii}     = f_loc;
            R2_tmp(ii)         = max(1 - mad(yf - Xf*B)/mad(yf), 0);
            dfe_tmp(ii)        = dfe_loc;
            Ball_tmp(:,ii)     = B + Priors_local;
        end

        % Scatter parfor results back into full-size arrays
        for ii = 1:nLst
            i = lst(ii);
            beta_all(:,i)     = beta_tmp(:,ii);
            tstat_all(:,i)    = tstat_tmp(:,ii);
            pval_all(:,i)     = pval_tmp(:,ii);
            ppos_all(:,i)     = ppos_tmp(:,ii);
            pneg_all(:,i)     = pneg_tmp(:,ii);
            P_all(i)          = P_tmp(ii);
            covb_all(:,:,i)   = covb_tmp(:,:,ii);
            w_all(:,i)        = w_tmp(:,ii);
            a_all{i}          = a_tmp{ii};
            sigma2_all(i)     = sigma2_tmp(ii);
            filter_all{i}     = filter_tmp{ii};
            R2_all(i)         = R2_tmp(ii);
            dfe_all(i)        = dfe_tmp(ii);
            Ball(:,i)         = Ball_tmp(:,ii);
        end

        values=unique(isRE(isRE~=0));
        for ii=1:length(values)
            Priors(isRE==values(ii)) = mean(Ball(isRE==values(ii),lst),2);

            Lambda(isRE==values(ii),isRE==values(ii))=diag(1./max(var(Ball(isRE==values(ii),lst),[],2),eps(1)));
            %Lambda(isRE==values(ii),isRE==values(ii))=1./cov(Ball(isRE==values(ii),:)');
            % otherwise, could use full covariance matrix here
        end

        outeriter=outeriter+1;

    end
end

% Assemble the stats struct from pre-allocated arrays
stats.beta   = beta_all;
stats.tstat  = tstat_all;
stats.dfe    = dfe_all(find(dfe_all, 1, 'first'));  % scalar (same for all channels)
stats.pval   = pval_all;
stats.ppos   = ppos_all;
stats.pneg   = pneg_all;
stats.P      = P_all;
stats.covb   = covb_all;
stats.w      = w_all;
stats.a      = a_all;
stats.sigma2 = sigma2_all;
stats.filter = filter_all;
stats.R2     = max(R2_all);

nirs.util.flushstdout(1);

end

%%
function out = myFilterLocal( f, y )
% here we are just making the first value zero before filtering to
% avoid weird effects introduced by zero padding
y1 = y(1,:);

y = bsxfun(@minus,y,y1);

out = filter(f, 1, y);
out = bsxfun(@plus,out,sum(f)*y1); % add the corrected offset back

end
