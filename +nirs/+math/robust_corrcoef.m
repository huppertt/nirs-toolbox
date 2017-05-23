function [r,p]=robust_corrcoef(d,verbose,mask);
% This is a robust version of the corrcoef function
% Based on Shevlyakov and Smirnov (2011). Robust Estimation of the
% Correlation Coeficient: An Attempt of Survey.  Austrian J. of Statistics.
%  40(1&2) pp. 147-156


if(nargin<2 || isempty(verbose))
    verbose=false;
end

if(nargin<3)
    mask=ones(size(d));
end


d=bsxfun(@minus,d,median(d,1));
d=bsxfun(@rdivide,d,sqrt(2)*mad(d,1,1));


W=ones(size(d,1),1);
for idx=1:size(d,2)
    W = min(W,wfun(d(:,idx)));
end
% Weight the whole model first to down-weight global motion artifacts
d=spdiags(W,0,length(W),length(W))*d;

% all pairwise combinations
r=ones(size(d,2));
n=size(d,2)^2;

if(verbose)
    disp('Progress');
    str='  0';
    fprintf('%s %%',str(end-2:end));
end


cnt=1;
for i=1:size(d,2)
    for j=1:size(d,2)
        if(verbose)
            str=['   ' num2str(round(100*cnt/n))];
            fprintf('\b\b\b\b\b%s %%',str(end-2:end));
        end
        warning('off','stats:statrobustfit:IterationLimit')
        lst=find(mask(:,i) & mask(:,j));
        
        d(lst,i)=bsxfun(@minus,d(lst,i),mean(d(lst,i)));
        d(lst,i)=bsxfun(@rdivide,d(lst,i),sqrt(var(d(lst,i))));
        
        d(lst,j)=bsxfun(@minus,d(lst,j),mean(d(lst,j)));
        d(lst,j)=bsxfun(@rdivide,d(lst,j),sqrt(var(d(lst,j))));
        r(i,j)=robustfit(d(lst,i),d(lst,j),'bisquare',[],'off');
        cnt=cnt+1;
    end
end
if(verbose)
    disp('completed');
end
% Section 2.3 Eqn 6
% r = sqrt(b1*b2)
r=sign(r).*abs(sqrt(r.*r'));
r=.5*(r+r');

n=sum(mask(:,1));

Tstat = r .* sqrt((n-2) ./ (1 - r.^2));
p = 2*tpvalue(-abs(Tstat),n-2);
p=tril(p,-1)+tril(p,-1)'+eye(size(p));

end

function w = wfun(r)
s = mad(r, 0) / 0.6745;
r = r/s/4.685;

w = (1 - r.^2) .* (r < 1 & r > -1);
end



function p = tpvalue(x,v)
%TPVALUE Compute p-value for t statistic.

normcutoff = 1e7;
if length(x)~=1 && length(v)==1
    v = repmat(v,size(x));
end

% Initialize P.
p = NaN(size(x));
nans = (isnan(x) | ~(0<v)); % v == NaN ==> (0<v) == false

% First compute F(-|x|).
%
% Cauchy distribution.  See Devroye pages 29 and 450.
cauchy = (v == 1);
p(cauchy) = .5 + atan(x(cauchy))/pi;

% Normal Approximation.
normal = (v > normcutoff);
p(normal) = 0.5 * erfc(-x(normal) ./ sqrt(2));

% See Abramowitz and Stegun, formulas 26.5.27 and 26.7.1.
gen = ~(cauchy | normal | nans);
p(gen) = betainc(v(gen) ./ (v(gen) + x(gen).^2), v(gen)/2, 0.5)/2;

% Adjust for x>0.  Right now p<0.5, so this is numerically safe.
reflect = gen & (x > 0);
p(reflect) = 1 - p(reflect);

% Make the result exact for the median.
p(x == 0 & ~nans) = 0.5;
end