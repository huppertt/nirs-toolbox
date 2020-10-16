function [r,p,W]=robust_corrcoef(d,verbose,mask);
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
mask = mask & ~isnan(d);
d(isnan(d)) = 0;

d=bsxfun(@minus,d,median(d,1));
d=bsxfun(@rdivide,d,1.4826*mad(d,1,1));


W=ones(size(d,1),1);
for idx=1:size(d,2)
    W = min(W,wfun(d(:,idx)));
end
% Weight the whole model first to down-weight global motion artifacts
d=spdiags(W,0,length(W),length(W))*d;

lst=find(mask(:,1));
[r,p]=corrcoef(d(lst,:));


p(:)=0;
m=size(r,1);
lowerhalf = (tril(ones(m),-1)>0);
rv = r(lowerhalf);
nv =sum(W);

% Tstat = +/-Inf and p = 0 if abs(r) == 1, NaN if r == NaN.
Tstat = rv .* sqrt((nv-2) ./ (1 - rv.^2));
p(lowerhalf) = 2*tpvalue(-abs(Tstat),nv-2);
p = p + p' + diag(diag(r)); % Preserve NaNs on diag.

end
% ------------------------------------------------
function w = wfun(r)
s = mad(r, 0) / 0.6745;
r = r/s/4.685;

w = (1 - r.^2) .* (r < 1 & r > -1);
end



% ------------------------------------------------
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
