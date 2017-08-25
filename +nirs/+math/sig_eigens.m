function lst=sig_eigens(s,pval)
% estimates the number of statistical eignevalues from an SVD decom

if(nargin<2)
    pval=0.05;
end

if(size(s,1)==size(s,2))
    s=diag(s);
end

err2 = 2*s.^2/length(s);

dL = s(1:end-1)-s(2:end);
t = dL ./sqrt(err2(1:end-1)+err2(2:end));

p=2*tcdf(-abs(t),length(s));

lst=1:max(find(p>pval));