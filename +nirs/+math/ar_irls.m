function stats = ar_irls( d,X,Pmax,tune )
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
    
    if nargin < 4
        tune = 4.685;
    end
    
    % preallocate stats
    nCond = size(X,2);
    nChan = size(d,2);
    nTime = size(d,1);
    
    stats.beta = zeros(nCond,nChan);    % betas
    stats.tstat = zeros(nCond,nChan);   % tstats
    stats.pval = zeros(nCond,nChan);    % two-sided t-test
    stats.ppos = zeros(nCond,nChan);    % one-sided t-test (positive only)
    stats.pneg = zeros(nCond,nChan);    % one-sided t-test (negative only)
    stats.P = zeros(nChan,1);           % the final AR model order
    stats.w = zeros(nTime,nChan);       % save the weights
    stats.dfe = nTime - nCond;          % degrees of freedom

%     
%        yfiltered=[];
%         weights=[];
%         Xfiltered=[];
    
    % loop through each channel
    for i = 1:nChan
        y = d(:,i);
        
        % initial fit
        B = X \ y;
        B0 = 1e6*ones(size(B));
        
        % iterative re-weighted least squares
        iter = 0;
        maxiter = 10;
        
        % while our coefficients are changing greater than some threshold
        % and it's less than the max number of iterations
        while norm(B-B0)/norm(B0) > 1e-2 && iter < maxiter
            % store the last fit
            B0 = B;
            
            % get the residual
            res = y - X*B;
                        
            % fit the residual to an ar model
            a = nirs.math.ar_fit(res, Pmax);
            
            % create a whitening filter from the coefficients
            f = [1; -a(2:end)];
            
            % filter the design matrix
            Xf = myFilter(f,X);
          
            % subtract constant from AR model and filter the data
            yf = myFilter(f,y);

            % perform IRLS
            [B, S] = robustfit(Xf,yf,'bisquare',tune,'off');
            
            iter = iter + 1;
        end
        
        % moco data & statistics
        stats.beta(:,i) = B;
        stats.tstat(:,i) = B./sqrt(diag(S.covb));
        stats.pval(:,i) = 2*tcdf(-abs(stats.tstat(:,i)),stats.dfe);     % two-sided
        stats.ppos(:,i) = tcdf(-stats.tstat(:,i),stats.dfe);            % one-sided (positive only)
        stats.pneg(:,i) = tcdf(stats.tstat(:,i),stats.dfe);             % one-sided (negative only)
        stats.P(i) = length(a)-1;
        
        L = pinv(Xf'*Xf); % more stable
        stats.covb(:,:,i) = L*S.mad_s^2;
        stats.w(:,i) = S.w;
        stats.a{i} = a;
        stats.sigma2(i)=S.mad_s^2;
        stats.filter{i}=f;
        stats.R2=max(1-mad(yf-Xf*B)/mad(yf),0);
%         yfiltered=[yfiltered; yf];
%         weights=[weights; S.w];
%         Xfiltered=sparse(blkdiag(Xfiltered,Xf));
        
    end   
    
   % diagnotics=fitglm(Xfiltered,yfiltered,'Weights',weights);
   % diagnotics=fitlm(Xfiltered,yfiltered,'Weights',weights);
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