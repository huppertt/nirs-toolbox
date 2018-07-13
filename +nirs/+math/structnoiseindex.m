function sni = structnoiseindex(d)
d=bsxfun(@minus,d,mean(d));
d=bsxfun(@minus,d,mean(d));

y=nirs.math.innovations(d,10);

sni=mad(d,1,1)./mad(y,1,1);