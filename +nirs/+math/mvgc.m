function [G, F, df1, df2, p] = mvgc(Y, Pmax, includeZeroLag, criterion)
% Multivariate granger causality
% [G, F, df1, df2, p] = mvgc(Y, Pmax , includeZeroLag );
if(nargin<4)
    criterion='AICc';
end
if(nargin<3)
    includeZeroLag=true;
end

if(includeZeroLag)
    s=0;
else
    s=1;
end
orders = s:Pmax;
q = length(orders);

% Center data
Y = bsxfun( @minus , Y , nanmean(Y) );

%% Create design matrix
[m, n, o] = size(Y); % [ time x channel x trial ]
X = zeros(m,n*q,o);

for i = 1:n   
    for j = 1:o
        inds = (i-1)*q+1 : i*q;
        X(:,inds,j) = nirs.math.lagmatrix( Y(:,i,j) , orders );
    end
end

lst = [reshape(repmat(1:n,[q 1]),[],1) repmat(orders',n,1)];

% Temporally concatenate trials
X = reshape( permute( X , [2 1 3] ) , [n*q m*o] )'; % [m n*q o] -> [m*o n*q]
Y = reshape( permute( Y , [2 1 3] ) , [n m*o] )';   % [m n o] -> [m*o n]
t = repmat( 1:m , 1 , o )';

if ~strcmpi(criterion,'MAX')
    %% Perform model order selection
    % Calculate log-likelihood for each model order
    % TODO - Behavior here currently doesn't reflect when includeZeroLag is true
    LogL = nan(Pmax,1);
    for i = 1:Pmax
        N = o*(m-i);
        time_inds = t>i;
        regr_inds = (lst(:,2)<=i);
        tmpX = X(time_inds,regr_inds);
        tmpY = Y(time_inds,:);
        b = tmpX \ tmpY;
        E = tmpY - tmpX*b;
        DSIG = det((E'*E)/(N-1));
        LogL(i) = -(N/2) * log(DSIG);
    end
    
    % Find model order with best information criterion
    crit = nirs.math.infocrit( LogL ,  o*(m-orders') , n^2 * orders' , criterion );
    Pmax = find(crit==nanmin(crit),1,'first');
    if isempty(Pmax), Pmax = 1; end

    % Prune higher model orders
    badlags = lst(:,2)>Pmax;
    X(:,badlags) = [];
    lst(badlags,:) = [];
end

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
G(G<0) = 0; % Due to numerical imprecision

df1 = Pmax * (n-1);          % Number of cross (non-auto) predictors [unrestricted parameters - restricted]
df2 = o*(m-Pmax) - (n*Pmax); % Effective observations (obs-order) - full # of predictors

assert(df2>0,'Degrees of freedom are too small. Increase number of samples or decrease model order.');

F = (exp(G)-1) .* (df2/df1); % F = (SSEr-SSEu)/df1 / SSEu/df2 = (SSEr-SSEu)/SSEu * df2/df1
                             % (SSEr-SSEu)/SSEu = SSEr/SSEu - 1 = exp(G) - 1
p = 1 - fcdf(F,df1,df2);

end


