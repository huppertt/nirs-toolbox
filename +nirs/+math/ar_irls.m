function [stats,resid] = ar_irls( d,X,Pmax,tune,nosearch,useGPU, singlePrecision)
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
    
    if nargin < 4 || isempty(tune)
        tune = 4.685;
    end
    if(nargin<5)
        nosearch=false;
    end
    if(nargin<6)
        useGPU=false;
    end

    if(nargin<7)
        singlePrecision=false;
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
        
        if ~any(y)  % In case y is all zeros (e.g., null hb)...
            %y = sqrt(0.01)*randn(length(y),1);  %...a non-null vector prevents a robustfit error without creating false positives
            Xfall{i}=nan(length(y),size(X,2));
            stats.beta(:,i) = nan(size(X,2),1);
            
            stats.covb(:,:,i) = nan(size(X,2),size(X,2));
            stats.w(:,i)=nan(length(y),1);
            stats.a{i} = [];
            stats.sigma2(i)=NaN;
            
            stats.tstat(:,i) =nan(size(X,2),1);
            stats.pval(:,i) = nan(size(X,2),1);     % two-sided
            stats.ppos(:,i) = nan(size(X,2),1);            % one-sided (positive only)
            stats.pneg(:,i) = nan(size(X,2),1);             % one-sided (negative only)
            
            resid(:,i)=nan(length(y),1);
            
            stats.filter{i}=[];
            stats.R2(i)=NaN;
        else
            
        % initial fit
        lstValid=~isnan(y);
        B = pinv(X(lstValid,:))* y(lstValid);
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
            if(length(Pmax)>1)
                p=Pmax(i);
            else
                p=Pmax;
            end
            a = nirs.math.ar_fit(res, p,nosearch);
            
            % create a whitening filter from the coefficients
            f = [1; -a(2:end)];
            
            % filter the design matrix
            Xf = myFilter(f,X);
          
            % subtract constant from AR model and filter the data
            lstInValid=isnan(y);
            lstValid=~isnan(y);
            if(nnz(lstInValid)>0)
                yy=y;
                yy(lstInValid)=interp1(find(lstValid),y(lstValid),find(lstInValid),'spline',true);
                 yf = myFilter(f,yy);
    
            else
            yf = myFilter(f,y);
            end
            % perform IRLS
            [B, S] = nirs.math.robustfit(Xf(lstValid,:),yf(lstValid),'bisquare',tune,'off');
            
            iter = iter + 1;
        end
        
        fprintf(1,'.');
        Xf=Xf(lstValid,:);
        
        % the model gets huge for EEG data.
        wXf=Xf;
        for id=1:size(wXf,2)
            wXf(:,id)=S.w.*wXf(:,id);
        end
        
        [U,~,~]=nirs.math.mysvd(wXf);
       %[U,~,~]=nirs.math.mysvd(diag(S.w)*Xf);
       
       % stats.dfe = length(yf)-sum(U(:).*U(:));
        stats.dfe = sum(S.w)-sum(U(:).*U(:));
        
        
        
        % note trace(A*B) = sum(reshape(A,[],1).*reshape(B',[],1)); 
        

        if(useGPU)
            
            %  Satterthwaite estimate of model DOF
            if(singlePrecision)
                g_Sw=gpuArray(single(S.w));
                g_wXf=gpuArray(single(wXf));
            else % double precision
                g_Sw=gpuArray(S.w);
                g_wXf=gpuArray(wXf);
            end

            gpuH=diag(g_Sw)-g_wXf*pinv(g_wXf'*g_wXf)*g_wXf';
            gpuHtH = gpuH' * gpuH;  

            % This order clears up memory using inplace insertion for
            % variables
            clear gpuH

            % Lower/denominator
            gpuHtH_sumsq=sum(reshape(gpuHtH,[],1).^2);
            clear gpuHtH
            
            % Upper/numerator
            gpuH=diag(g_Sw)-g_wXf*pinv(g_wXf'*g_wXf)*g_wXf';
            gpuH_upper=sum(reshape(gpuH,[],1).*reshape(gpuH',[],1))^2;
            clear gpuH

            stats.dfe =gather(gpuH_upper/gpuHtH_sumsq);

            %stats.dfe =gather(sum(reshape(gpuH,[],1).*reshape(gpuH',[],1))^2/sum(reshape(gpuHtH,[],1).^2));
            
        else
            %  Satterthwaite estimate of model DOF
            if(singlePrecision)
                sSw=single(S.w);
                swXF=single(wXf);
                H=diag(sSw)-swXF*pinv(swXF'*swXF)*swXF';
                HtH=H'*H;
            else % double precision
                H=diag(S.w)-wXf*pinv(wXf'*wXf)*wXf';
                HtH=H'*H;
            end
            
            stats.dfe =sum(reshape(H,[],1).*reshape(H',[],1))^2/sum(reshape(HtH,[],1).^2);

        end

        
        
        %stats.dfe = trace(H'*H)^2/trace(H'*H*H*H');  % same result but 4-6x slower
        
        % moco data & statistics
        stats.beta(:,i) = B;
        stats.P(i) = length(a)-1;
        
        L = pinv(Xf'*Xf); % more stable
        Xfall{i}=nan(length(y),size(wXf,2));
        Xfall{i}(lstValid,:)=wXf;
        stats.covb(:,:,i) = L*S.robust_s^2;  
        stats.w(:,i)=nan(length(y),1);
        stats.w(lstValid,i) = S.w;
        stats.a{i} = a;
        stats.sigma2(i)=S.robust_s^2;  
        
        stats.tstat(:,i) = stats.beta(:,i)./sqrt(diag(S.covb));
        stats.pval(:,i) = 2*tcdf(-abs(stats.tstat(:,i)),stats.dfe);     % two-sided
        stats.ppos(:,i) = tcdf(-stats.tstat(:,i),stats.dfe);            % one-sided (positive only)
        stats.pneg(:,i) = tcdf(stats.tstat(:,i),stats.dfe);             % one-sided (negative only)

        resid(:,i)=nan(length(y),1);
        resid(lstValid,i)=S.resid.*S.w; 
        
        stats.filter{i}=f;
        sse =  norm(yf(lstValid) - Xf*B)^2;
         
        sst =  norm(yf(lstValid) - mean(yf(lstValid)))^2;
        stats.R2(i)=1-sse./sst;
        
        N=length(y)-size(B,2)-stats.P(i);
        stats.logLik(i)=(-N/2)*log(2*pi*sse/N)-N/2;
       
        
        end
        
    end   
   
     
    resid=resid-ones(size(resid,1),1)*nanmedian(resid,1);
    if(size(resid,2)>400)
        C=resid'*resid;
    else
        for i=1:size(stats.beta,2)
            for j=i:size(stats.beta,2)
                a=resid(:,i);
                b=resid(:,j);
                C(i,j)=1.4810*nanmedian(a.*b);  % var(x) = 1.4810 * MAD(x,1)
                C(j,i)=C(i,j);
            end

        end
    end
    C=C*(nanmean(stats.sigma2'./diag(C)));   %fix the scaling due to the dof (which is a bit hard to track because it changes per channel, so use the average)
    
    
    if(length(Xfall)>200)
        covb=zeros(size(stats.beta,1),size(stats.beta,1),size(stats.beta,2),size(stats.beta,2));
    
         Xf=zeros(size(Xfall{1})); for i=1:length(Xfall); Xf=Xf+Xfall{i}; end;
         Xf=Xf/length(Xfall);
         
         ixtx=pinv(Xf'*Xf);
         for i=1:size(stats.beta,2)
           for j=i:size(stats.beta,2)
               covb(:,:,i,j) =covb(:,:,i,j)+ixtx*C(i,j);
               covb(:,:,j,i) =covb(:,:,j,i)+ixtx'*C(j,i);  % done to ensure symmetry
            end
        end

    else
        covb=zeros(size(stats.beta,1),size(stats.beta,1),size(stats.beta,2),size(stats.beta,2));
        for i=1:size(stats.beta,2)
            for j=1:size(stats.beta,2)
                lstV=~isnan(sum(Xfall{i},2)+sum(Xfall{j},2));
                ixtx=pinv(Xfall{i}(lstV,:)'*Xfall{j}(lstV,:));
                covb(:,:,i,j) =covb(:,:,i,j)+ixtx*C(i,j);
                covb(:,:,j,i) =covb(:,:,j,i)+ixtx'*C(j,i);  % done to ensure symmetry
            end
        end
    end
    covb=covb/2;

    
%     
%     figure(1); cla; d=[];
%     for i=1:32; d(i)=squeeze(covb(2,2,i,i)); end;
%     plot(log(squeeze(stats.covb(2,2,:))),'b')
%     hold on;
%     plot(log(d),'r--')
%     pause;
%     
    stats.covb = real(covb);
% 
%     for i=1:size(stats.beta,2)
%         stats.tstat(:,i) = stats.beta(:,i)./sqrt(diag(squeeze(stats.covb(:,:,i,i))));
%         stats.pval(:,i) = 2*tcdf(-abs(stats.tstat(:,i)),stats.dfe);     % two-sided
%         stats.ppos(:,i) = tcdf(-stats.tstat(:,i),stats.dfe);            % one-sided (positive only)
%         stats.pneg(:,i) = tcdf(stats.tstat(:,i),stats.dfe);             % one-sided (negative only)
%     end
    
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