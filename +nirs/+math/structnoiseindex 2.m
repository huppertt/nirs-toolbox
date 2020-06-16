function [sni,moco] = structnoiseindex(d,n)
d=bsxfun(@minus,d,mean(d));
d=bsxfun(@minus,d,mean(d));

if(nargin<2)
    n=5;
end

y=nirs.math.innovations(d,n);

tiny_s = 1e-6 * std(y);
t = 4.685;
s = median(abs(y)) / 0.6745;
r=y./(ones(size(y,1),1)*(max(s,tiny_s)*t));
w = (abs(r)<1) .* (1 - r.^2).^2;
w=sqrt(w);
   
sni=mad(d,1,1)./mad(y,1,1);

moco=sum(w,1)/length(w);