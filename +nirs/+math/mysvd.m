function [u,s,v]=mysvd(X)

if(isdiag(X) & issparse(X))
    u=speye(size(X));
    v=speye(size(X));
    s=X;
    return;
end

if(issparse(X)); X=full(X); end;

if(size(X,1)>size(X,2))
    [u,s,v]=svd(X,0);
else
    [v,s,u]=svd(X',0);
end

return