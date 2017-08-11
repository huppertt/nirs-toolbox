function [r,p] = robust_xcorrcoef( d , maxlags , verbose , mask )

if(nargin<2 || isempty(maxlags))
    maxlags=size(d,1)/2;
end

if(nargin<3 || isempty(verbose))
    verbose=false;
end

if(nargin<4)
    mask=ones(size(d));
end
mask = mask & ~isnan(d);
d(isnan(d)) = 0;

[M,N] = size(d);
maxlags = ceil(maxlags);
lags = -maxlags:maxlags;
numlag = length(lags);
r = zeros(N,N,numlag);

d=bsxfun(@minus,d,median(d,1));
d=bsxfun(@rdivide,d,1.4826*mad(d,1,1));

W=ones(size(d,1),1);
for idx=1:size(d,2)
    W = min(W,wfun(d(:,idx)));
end
% Weight the whole model first to down-weight global motion artifacts
d=spdiags(W,0,length(W),length(W))*d;

warning('off','stats:statrobustfit:IterationLimit');

for i = 1:N
    for j = i:N
        
        lst=(mask(:,i) & mask(:,j));
        tmp_d = d(lst,[i j]);
        tmp_M = size(tmp_d,1);
        tmp_inds = 1:tmp_M;
        
        for k = 1:numlag
                        
            lag = lags(k);
            iinds = tmp_inds;
            oinds = tmp_inds - lag;
            dropinds = oinds<1 | oinds>tmp_M;
            iinds(dropinds)=[];
            oinds(dropinds)=[];
            
            d1 = tmp_d(:,1);
            d2 = zeros(tmp_M,1);
            d2(oinds) = tmp_d(iinds,2);
            
            b1 = robustfit(d1,d2,'bisquare',[],'on');
            b2 = robustfit(d2,d1,'bisquare',[],'on');
            
            val = sign(b1(2)+b2(2)) * sqrt(abs(b1(2).*b2(2)));
            
            r(i,j,numlag-k+1) = val;
            r(j,i,k) = val;
        end
    end
end

n=sum(mask(:,1));

Tstat = r .* sqrt((n-2) ./ (1 - r.^2));
p = 2*nirs.math.tpvalue(-abs(Tstat),n-2);
p=p+repmat(eye(N),[1 1 numlag]);

end

function w = wfun(r)
s = mad(r, 0) / 0.6745;
r = r/s/4.685;

w = (1 - r.^2) .* (r < 1 & r > -1);
end
