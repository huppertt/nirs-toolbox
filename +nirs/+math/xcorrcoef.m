function [r,p]=xcorrcoef(d,maxlags,mask)
% This returns the largest magnitude cross-correlation coefficient
% [r,p]=nirs.math.xcorrcoef(d,maxlags,mask)

if(nargin<2 || isempty(maxlags))
    maxlags=size(d,1);
end

if(ischar(maxlags))
    maxlags = Fs*str2double(maxlags(1:strfind(maxlags,'x')-1));
end

if(nargin<3)
    mask=ones(size(d));
end

% all pairwise combinations
r=ones(size(d,2));
n=size(d,2)^2;

cnt=1;
for i=1:size(d,2)
    for j=1:size(d,2)
        lst=find(mask(:,i) & mask(:,j));
        rs = xcorr(d(lst,i),d(lst,j),maxlags,'Coeff');
        r(i,j)=rs(find(abs(rs)==max(abs(rs)),1));
        cnt=cnt+1;
    end
end
r(r>1) = 1;
r(r<-1) = -1;

n=size(d,1);

Tstat = r .* sqrt((n-2) ./ (1 - r.^2));
p = 2*tpvalue(-abs(Tstat),n-2);
p=tril(p,-1)+tril(p,-1)'+eye(size(p));

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