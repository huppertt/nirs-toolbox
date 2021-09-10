function [stats,resid] = ar_irls_ftest(y,X,p,tune)

    if nargin < 4 || isempty(tune)
        tune = 4.685;
    end

    [stats,resid]=nirs.math.ar_irls(y,X,p,tune);
    for i=1:size(X,2)
        lst=1:size(X,2);
        lst(i)=[];
        Stats2=nirs.math.ar_irls(y,X(:,lst),stats.P,tune,true);
        LL=-2*[Stats2.logLik'-stats.logLik'];
        stats.Fpval(i,:)=1-chi2cdf(LL,1);
    end
    
  

end
