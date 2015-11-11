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

for j = 1:n
    for i = 1:n
        
        % restricted model
        rLst = find(lst(:,1)~=j & ~(lst(:,1)==i & lst(:,2)==0));
        %unrestricted model
        uLst = find( ~(lst(:,1)==i & lst(:,2)==0));
        
        iXu = inv(X(:,uLst)'*X(:,uLst))*X(:,uLst)';
        iXr = inv(X(:,rLst)'*X(:,rLst))*X(:,rLst)';
        
        % unrestricted model (all terms)
        a   = iXu * Y(:,i);
        ru  = Y(:,i) - X(:,uLst)*a;
        
        % restricted model (no j terms)
        a  = iXr * Y(:,i);
        rr = Y(:,i) - X(:,rLst)*a;
        
        G(i,j)      = log(mad(rr)/mad(ru));
        df1(i,j)    = size(X,2)-length(rLst);
        df2(i,j) 	= size(X,1) - size(X,2);
        F(i,j)      = ((mad(rr)^2-mad(ru)^2)/df1(i,j))/(mad(ru)^2/df2(i,j));
        F(i,j)      = max(F(i,j),0);
        p(i,j)      = fcdf(1/F(i,j), df2(i,j), df1(i,j));
        
    end
    
end


