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

% model selection
Xall = []; Yall = [];
for i = 1:n
    chan_inds = lst(:,1) == i;
    Xall = [Xall; X(:,chan_inds)];
    Yall = [Yall; Y(:,i)];
end
[~, ~, crit] = nirs.math.stepwise( Xall , Yall , 'BIC' );
Pmax = find(crit==min(crit),1,'first');
badlags = lst(:,2)>Pmax;
X(:,badlags) = [];
lst(badlags,:) = [];

G = nan(n);

% Create unrestricted models (same for all channels if no zero lag)
if ~includeZeroLag
    iXu = inv(X'*X)*X';
end

% Create restricted models (1 per channel if no 0-lag, 1/channel^2 if 0-lag)
% Todo: A better way for 0-lag is to compute the inv(Xr'*Xr) of restricted model once for each channel,
%       then remove the effect of self 0-lag as needed (https://emtiyaz.github.io/Writings/OneColInv.pdf)
if ~includeZeroLag
    iXrs = zeros([size(X,2)-Pmax size(X,1) n]);
    for j = 1:n
        % restricted model
        LstR = lst(:,1)~=j;
        Xr = X(:,LstR);

        iXrs(:,:,j) = inv(Xr'*Xr)*Xr'; 
    end
end

for i = 1:n
    
    % solve unrestricted model (all terms)
    if includeZeroLag
        
        uind = ( ~(lst(:,1)==i & lst(:,2)==0));
        LstU = lst(uind,:);
        Xu = X(:,uind);

        iXu = inv(Xu'*Xu)*Xu';
        
    else
        LstU = lst;
        Xu = X;
    end
    
    a   = iXu * Y(:,i);
    ru  = Y(:,i) - Xu*a;
    
    SSEu = var(ru);
    
    for j = 1:n
    
        % restricted model
        LstR = (LstU(:,1)~=j);
        Xr = Xu(:,LstR);
        
        if includeZeroLag
            iXr = inv(Xr'*Xr)*Xr';
        else
            iXr = iXrs(:,:,j);
        end
        
        % restricted model (no j terms)
        a  = iXr * Y(:,i);
        rr = Y(:,i) - Xr*a;
        
        SSEr = var(rr);
        
        G(i,j) = log(SSEr/SSEu);

    end
    
end

df1 = Pmax * (n-1);          % Number of cross (non-auto) predictors [unrestricted parameters - restricted]
df2 = (m-Pmax) - (n*Pmax+1); % Effective observations (obs-order) - full # of predictors

F = (exp(G)-1) .* (df2/df1); % F = (SSEr-SSEu)/df1 / SSEu/df2 = (SSEr-SSEu)/SSEu * df2/df1
                             % (SSEr-SSEu)/SSEu = SSEr/SSEu - 1 = exp(G) - 1
p = 1 - fcdf(F,df1,df2);

end


