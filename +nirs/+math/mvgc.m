function [G, F, df1, df2, p] = mvgc(Y, Pmax, includeZeroLag)
% Multivariate granger causality
% [G, F, df1, df2, p] = mvgc(Y, Pmax , includeZeroLag );
 
if(nargin<3)
    includeZeroLag=true;
end

if(includeZeroLag)
    s=0;
else
    s=1;
end

%% Create design matrix
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

X = [ones(size(X,1),1) X];
lst = [0 0; lst];

%% Perform model order selection
% Calculate log-likelihood for each model order
% TODO - Behavior here currently doesn't reflect when includeZeroLag is true
LogL = nan(Pmax,1);
for i = 1:Pmax
    time_inds = (i+1):m;
    regr_inds = (lst(:,2)<=i);
    tmpX = X(time_inds,regr_inds);
    tmpY = Y(time_inds,:);
    b = tmpX \ tmpY;
    E = tmpY - tmpX*b;
    DSIG = det((E'*E)/(m-i-1));
    LogL(i) = -((m-i)/2) * log(DSIG);
end

% Find model order with best information criterion
crit = nirs.math.infocrit( LogL ,  m-(1:Pmax)' , n^2 * (1:Pmax)' , 'AICc' );
Pmax = find(crit==nanmin(crit),1,'first');

% Prune higher model orders
badlags = lst(:,2)>Pmax;
X(:,badlags) = [];
lst(badlags,:) = [];

% Remove early timepoints that don't have predictive information to model
X(1:Pmax,:) = [];
Y(1:Pmax,:) = [];

%% Calculate restricted and unrestricted OLS model fits
if ~includeZeroLag
    
    % Unrestricted model
    b = X \ Y;
    E = Y - X*b;
    SSEu = var(E)';
    
    % Restricted model
    SSEr = nan(n);
    for i=1:n
        LstR = lst(:,1)~=i;
        b = X(:,LstR) \ Y;
        E = Y - X(:,LstR)*b;
        SSEr(:,i) = var(E)';
    end
    
else
    
    SSEu = nan(n,1);
    SSEr = nan(n);
    
    for i = 1:n
        
        % Unrestricted model
        LstU = ~(lst(:,1)==i & lst(:,2)==0);
        Xu = X(:,LstU);

        b = Xu \ Y(:,i);
        E = Y(:,i) - Xu*b;
        SSEu(i) = var(E)';
        
        % Restricted model
        for j = 1:n
            
            LstR = lst(:,1)~=j & ~(lst(:,1)==i & lst(:,2)==0);
            Xr = X(:,LstR);

            b = Xr \ Y(:,i);
            E = Y(:,i) - Xr*b;
            SSEr(i,j) = var(E)';
            
        end
    end
    
end

%% Compute granger stats
G = bsxfun( @minus , log(SSEr) , log(SSEu) );

df1 = Pmax * (n-1);          % Number of cross (non-auto) predictors [unrestricted parameters - restricted]
df2 = (m-Pmax) - (n*Pmax);   % Effective observations (obs-order) - full # of predictors

F = (exp(G)-1) .* (df2/df1); % F = (SSEr-SSEu)/df1 / SSEu/df2 = (SSEr-SSEu)/SSEu * df2/df1
                             % (SSEr-SSEu)/SSEu = SSEr/SSEu - 1 = exp(G) - 1
p = 1 - fcdf(F,df1,df2);

end


