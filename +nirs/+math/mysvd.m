function [u,s,v]=mysvd(X)

if(size(X,1)>size(X,2))
    [u,s,v]=svd(X,0);
else
    [v,s,u]=svd(X',0);
end


return