function  [Bmat,Cov] = linsolveW(Xmat,Ymat,Cov,pmax,opt);
% this does a weighted, robust, AR-whitened version of linsolve

W = inv(cholcov(Cov));

XW = W*Xmat;
YW = W*Ymat;

if(pmax==0)
    Bmat=linsolve(XW,YW);
    Cov=cov(YW-XW*Bmat);
else
    stats = nirs.math.ar_irls(YW,XW,pmax);
    Bmat=stats.beta;
    Cov=stats.covb;
end



return
