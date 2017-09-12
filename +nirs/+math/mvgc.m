function [G, F, df1, df2, p] = mvgc(Y, Pmax,includeZeroLag)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
 
if(nargin<3)
    includeZeroLag=true;
end

Y=Y-(ones(size(Y,1),1)*mean(Y,1));
Y=Y-(ones(size(Y,1),1)*mean(Y,1));


 if(includeZeroLag)
     s=0;
 else
     s=1;
 end
 
[m, n] = size(Y);
X = []; lst = [];
% model selection
for i = 1:n
    
    % column indices for channel i
    ll=s:Pmax;
    idx = size(lst,1)+1:size(lst,1)+length(ll);    
    % keep track of them
    lst(idx,1) = i;
    lst(idx,2) = ll;
    
    % design matrix
    X(:,idx) = nirs.math.lagmatrix(Y(:,i), s:Pmax);
end




%Remove an DC shifts in the data
l=find(~all(X==ones(size(X,1),1)*X(1,:)));
X=X(:,l);

% add a constant term
X   = [ones(size(X,1),1) X];
lst = [0 0; lst];

[G,F,df1,df2,p] = deal(nan(n));

for i = 1:n

    % unrestricted model
    uind = ( ~(lst(:,1)==i & lst(:,2)==0));
    LstU = lst(uind,:);
    Xu = X(:,uind);
    
    iXu = inv(Xu'*Xu)*Xu';
    
    % unrestricted model (all terms)
    a   = iXu * Y(:,i);
    ru  = Y(:,i) - Xu*a;
    
    SSEu = mad(ru)^2;
    
    for j = 1:n
    
        % restricted model
        LstR = (LstU(:,1)~=j);
        Xr = Xu(:,LstR);
        
        iXr = inv(Xr'*Xr)*Xr';
        
        % restricted model (no j terms)
        a  = iXr * Y(:,i);
        rr = Y(:,i) - Xr*a;
        
        SSEr = mad(rr)^2;
        
        G(i,j)      = log(SSEr/SSEu);
        df1(i,j)    = size(Xu,2)-sum(LstR);
        df2(i,j)    = size(Xu,1)-length(ll)-size(Xu,2);
        F(i,j)      = ((SSEr-SSEu)/df1(i,j))/(SSEu/df2(i,j));
        F(i,j)      = max(F(i,j),0);
        p(i,j)      = fcdf(1/F(i,j), df2(i,j), df1(i,j));
        
    end
    
end

end


