function [G, F, df1, df2, p] = roubst_mvgc(Y, Pmax,includeZeroLag)
% Robust version of multi-variate Grangers

   niter=1;
   

%Do it twice to avoid numerical issues
Y = Y-ones(size(Y,1),1)*mean(Y,1);
Y = Y-ones(size(Y,1),1)*mean(Y,1);


w = ones(size(Y,1),1);
for idx=1:size(Y,2)
    w=min(w,wfun(Y(:,idx)));
end
W=diag(w);



%%

if(nargin<3)
    includeZeroLag=true;
end


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

Y= W*Y;
X= W*X;



for j = 1:n
    for i = 1:n
        
        % restricted model
        rLst = find(lst(:,1)~=j & ~(lst(:,1)==i & lst(:,2)==0));
        %unrestricted model
        uLst = find( ~(lst(:,1)==i & lst(:,2)==0));
        
        
       
        % % unrestricted model (all terms)
        iXu = inv(X(:,uLst)'*X(:,uLst))*X(:,uLst)';
        a   = iXu * Y(:,i);
           
      
        for iter=1:niter
            w = diag(wfun(Y(:,i)-X(:,uLst)*a));
            wX=w*X(:,uLst);
            iXu = inv(wX'*wX)*wX';
            a   = iXu * Y(:,i);
        end
        %a=robustfit(X(:,uLst),Y(:,i),'bisquare',4.685,'off');
        ru  = Y(:,i) - X(:,uLst)*a;
        
        % % restricted model (no j terms)
         iXr = inv(X(:,rLst)'*X(:,rLst))*X(:,rLst)';
         a  = iXr * Y(:,i);
        
      
        for iter=1:niter
            w = diag(wfun(Y(:,i)-X(:,rLst)*a));
            wX=w*X(:,rLst);
            iXr = inv(wX'*wX)*wX';
            a   = iXr * Y(:,i);
        end
       % a=robustfit(X(:,rLst),Y(:,i),'bisquare',4.685,'off');
       rr = Y(:,i) - X(:,rLst)*a;
        
        G(i,j)      = log(mad(rr)/mad(ru));
        df1(i,j)    = size(X,2)-length(rLst);
        df2(i,j) 	= size(X,1) - size(X,2);
        F(i,j)      = ((mad(rr)^2-mad(ru)^2)/df1(i,j))/(mad(ru)^2/df2(i,j));
        F(i,j)      = max(F(i,j),0);
        p(i,j)      = fcdf(1/F(i,j), df2(i,j), df1(i,j));
        
    end
    
    
end

end



function w = wfun(r)
    s = mad(r, 0) / 0.6745;
    r = r/s/4.685;
    
    w = (1 - r.^2) .* (r < 1 & r > -1);
end