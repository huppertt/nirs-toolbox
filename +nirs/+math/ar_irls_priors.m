function stats = ar_irls_priors( d,X,isRE,Pmax,tune,type )
% isFE=[1 0]; -  boolean array of if Fixed Effects (or random effects) in
% the model- defines which regressors have spatial covariance priors
% applied to them

warning('off','stats:statrobustfit:IterationLimit')

if nargin < 5
    tune = 4.685;
end

Priors=zeros(size(X,2),1);
Lambda=zeros(size(X,2));

nChan=size(d,2);



utype=unique(type);
ntypes=length(utype);

for it=1:ntypes
    lst=find(ismember(type,utype{it}));
    
    Lambda0=1E4*eye(size(X,2));
    outeriter=0;
    
    while(norm(Lambda0-Lambda)>1E-4 & outeriter<5)
        Lambda0=Lambda;
        for ii = 1:length(lst)
            i=lst(ii);
            disp([num2str(i) ' of ' num2str(nChan)]);
            y = d(:,i)-X*Priors;
            B = pinv(X'*X+Lambda)*X'*y;
            
            for iter2=1:5
                res = y - X*B;
                a = nirs.math.ar_fit(res, Pmax);
                f = [1; -a(2:end)];
                
                % filter the design matrix
                Xf = myFilter(f,X);
                
                % subtract constant from AR model and filter the data
                yf = myFilter(f,y);
                
                % [B, S] = robustfit(Xf,yf,'bisquare',tune,'off');
                for iter=1:10
                    r = yf - Xf*B;
                    s = mad(r, 0) / 0.6745;
                    r = r/s/tune;
                    w = (1 - r.^2) .* (r < 1 & r > -1);
                    Xfw=diag(w)*Xf;
                    yfw=diag(w)*yf;
                    B = pinv(Xfw'*Xfw+Lambda)*Xfw'*yfw;
                end
                
                
                r = yf - Xf*B;
                covb=pinv(Xfw'*Xfw+Lambda)*mean(r.^2);
                stats.beta(:,i) = B+Priors;
                stats.tstat(:,i) = B./sqrt(diag(covb));
                stats.tstat=real(stats.tstat);
                stats.tstat(isnan(stats.tstat))=0;
                stats.dfe=length(yf)-length(B);
                stats.pval(:,i) = 2*tcdf(-abs(stats.tstat(:,i)),stats.dfe);     % two-sided
                stats.ppos(:,i) = tcdf(-stats.tstat(:,i),stats.dfe);            % one-sided (positive only)
                stats.pneg(:,i) = tcdf(stats.tstat(:,i),stats.dfe);             % one-sided (negative only)
                stats.P(i) = length(a)-1;
                
                stats.covb(:,:,i) = covb;
                stats.w(:,i) = w;
                stats.a{i} = a;
                stats.sigma2(i)=mad(r)^2;
                stats.filter{i}=f;
                stats.R2=max(1-mad(yf-Xf*B)/mad(yf),0);
                
                Ball(:,i)=B+Priors;
            end
            
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