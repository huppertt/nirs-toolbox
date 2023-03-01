function [stats,resid] = ar_irls_ftest(y,X,p,tune,useGPU,singlePrecision)

    if nargin < 4 || isempty(tune)
        tune = 4.685;
    end

    if nargin <5
        useGPU=false;
    end

    if nargin <6
        singlePrecision=false;
    end

    [stats,resid]=nirs.math.ar_irls(y,X,p,tune,false,useGPU,singlePrecision);
    for i=1:size(X,2)
        lst=1:size(X,2);
        lst(i)=[];
        Stats2=nirs.math.ar_irls(y,X(:,lst),stats.P,tune,false,useGPU,singlePrecision);
        LL=-2*[Stats2.logLik'-stats.logLik'];
        stats.Fpval(i,:)=1-chi2cdf(LL,1);
    end
    
  

end
