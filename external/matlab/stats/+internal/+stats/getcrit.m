function [crit,pval] = getcrit(ctype, alpha, df, ng, t)
%GETCRIT Compute critical value for multiple comparisons.
%        Also return the p-value for the t statistic.

%   Copyright 1993-2013 The MathWorks, Inc.

crit = Inf;
if nargin<5
    t = [];
end

kwds = {'tukey-kramer' 'dunn-sidak' 'bonferroni' 'scheffe' 'lsd' 'hsd'};
ctypes = internal.stats.getParamVal(ctype,kwds,'''ComparisonType''',true);

for j=1:length(ctypes)
    onetype = ctypes{j};
   
   switch onetype
    case {'tukey-kramer' 'hsd'}
     crit1 = internal.stats.stdrinv(1-alpha, df, ng) / sqrt(2);
     
     % The T-K algorithm is inaccurate for small alpha, so compute
     % an upper bound for it and make sure it's in range.
     ub = internal.stats.getcrit('dunn-sidak', alpha, df, ng);
     if (crit1 > ub), crit1 = ub; end
     pval = 0*t;
     for k=1:numel(t)
         pval(k) = 1 - internal.stats.stdrcdf(sqrt(2)*abs(t(k)),df,ng);
     end

    case 'dunn-sidak'
     kstar = nchoosek(ng, 2);
     alf = 1-(1-alpha).^(1/kstar);
     crit1 = tinv(1-alf/2, df);
     pval = 1 - (1-2*tcdf(-abs(t),df)).^kstar;

    case 'bonferroni'
     kstar = nchoosek(ng, 2);
     crit1 = -tinv(alpha / (2*kstar), df);
     pval = 2*kstar*tcdf(-abs(t),df);

    case 'lsd'
     crit1 = -tinv(alpha / 2, df);
     pval = 2*tcdf(-abs(t),df);

    case 'scheffe'
     tmp = finv(1-alpha, ng-1, df);
     crit1 = sqrt((ng-1) * tmp);
     pval = fcdf((t.^2)/(ng-1),ng-1,df,'upper');
     
    otherwise
     error(message('stats:multcompare:BadCType', ctype));
   end

   if (~isnan(crit1))
       crit = min(crit, crit1);
   end
end
pval = min(1,pval); % multiple comparison may push out of range
