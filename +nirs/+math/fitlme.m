function [beta,bHat,covb,LL,w] = fitlme(X,Y,Z,robust_flag,zero_theta,verbose)
% Robust linear mixed-effects model fiting
% [beta,bHat,covb,LL,w] = nirs.math.fitlme( X, Y, Z, robust_flag, zero_theta_flag, verbose_flag )
%
% TODO: Add support for covariance patterns other than isotropic
if nargin<4, robust_flag = false; end
if nargin<5, zero_theta = false; end
if nargin<6, verbose = false; end

nT = size(Y,1);
nY = size(Y,2);
nX = size(X,2);
nZ = size(Z,2);

%% Handle all-NaN variables
bad_vars = all(isnan(Y),1) | all(bsxfun(@eq,Y,Y(1,:)),1);
if any(bad_vars)
    
    beta = nan(nX,nY);
    bHat = nan(nZ,nY);
    covb = nan(nX,nX,nY);
    LL = nan(1,nY);
    w = nan(nT,nY);
    
    [beta(:,~bad_vars),bHat(:,~bad_vars),covb(:,:,~bad_vars),LL(~bad_vars),w(:,~bad_vars)] = nirs.math.fitlme(X,Y(:,~bad_vars),Z,robust_flag,zero_theta,verbose);
    
    return;
end

%% Handle bad time points
bad_times = any(isnan(X),2) | any(isnan(Y),2) | any(isnan(Z),2);
if any(bad_times)
    X(bad_times,:) = [];
    Y(bad_times,:) = [];
    Z(bad_times,:) = [];
end
nT0 = nT;
nT = size(X,1);
beta = nan(nX,nY);
bHat = nan(nZ,nY);
covb = nan(nX,nX,nY);
LL = nan(1,nY);
w = nan(nT,nY);

if isempty(X) || isempty(Y), return; end

%% Run separately on each unique predictee
if nY > 1
    [~,uinds,indsu] = unique(Y','rows','stable');

    for i = 1:length(uinds)

        ind = uinds(i);
        out = find(indsu == i);
        [ubeta,ubHat,ucovb,uLL,uw] = nirs.math.fitlme(X,Y(:,ind),Z,robust_flag,zero_theta,verbose);
        
        beta(:,out) = repmat(ubeta,[1 length(out)]);
        bHat(:,out) = repmat(ubHat,[1 length(out)]);
        covb(:,:,out) = repmat(ucovb,[1 1 length(out)]);
        LL(out) = uLL;
        w(~bad_times,out) = repmat(uw,[1 length(out)]);
        
    end
    
    return;
end

%% Ensure no sparse matrices
X = full(X);
Y = full(Y);
Z = full(Z);

%% Compute theta
if zero_theta
    theta = 0;
else
    theta = solveForTheta(X,Y,Z);
end

%% Solve initial model
[beta,bHat,covb,LL,sigma2] = solveLME(X,Y,Z,theta);

%% Robust loop
if robust_flag

    iter = 1;
    tune = 4.685;
    D = sqrt(eps(class(X)));

    % Adjust by leverage to account for prior weight differences in design matrix
    lev = diag( full(X) * pinv(full(X)) );
    adj = 1 ./ sqrt(1-min(.9999,lev));
    xrank = rank(X);
    num_params = max(1,xrank);
    
    while iter<50

        % Calculate residuals and weights from previous iteration
        resid = (Y - X*beta - Z*bHat) .* adj;
        resid_s = studentizeResiduals( resid , num_params );
        w = bisquare( resid_s , tune );
        w = diag(w);
        beta0=beta;
        
        % Bail out if weights are bad
        if length(unique(w*Y))<2 || any(isnan(w(:))) || any(~isfinite(w(:)))
            beta(:)=nan;
            bHat(:)=nan;
            covb(:)=nan;
            LL(:)=nan;
            return;
        end
        
        % Re-estimate using new weights
        if ~zero_theta
            theta = solveForTheta(w*X,w*Y,w*Z,theta); % Get optimal theta
        end
        [beta,bHat] = solveLME(w*X,w*Y,w*Z,theta); % Solve model

        if verbose
            disp(['Robust fit iteration ' num2str(iter) ' : ' num2str(max(abs(beta-beta0)))]);
        end
        
        % Terminate if estimated coefficients have converged
        if ~any(abs(beta-beta0) > D*max(abs(beta),abs(beta0)))
            break
        end
        
        iter=iter+1;
    end
    
    % Calculate sigma
    sigma_robust = robustSigma(X,Y,Z,adj,num_params,tune,beta,bHat);
    sigma = max( sigma_robust , sqrt((sigma2 * xrank^2 + sigma_robust^2 * nT) / (xrank^2 + nT)) );
    
    % Calculate covariance betas
    Lambda = sqrt(exp(theta)) * speye(nZ); % Isotropic covariance pattern
    Iq = spdiags(ones(nZ,1),0,nZ,nZ);
    [R,~,S] = chol( Lambda'*sparse(Z'*Z)*Lambda + Iq );
    Q1 = ((X'*Z*Lambda)*S) / R;
    R1R1t = X'*X - Q1*Q1';
    R1 = cholSafe(R1R1t,'lower');
    invR1 = R1(1:xrank,1:xrank) \ eye(xrank);
    covb = sigma^2*(invR1'*invR1);
    w0 = diag(w);
    w = zeros(nT0,1);
    w(~bad_times) = w0;
else
    w = ones(nT0,1);
end

end

%% Use unconstrained nonlinear optimization to find theta that minimizes log-likelihood
function theta = solveForTheta(X,y,Z,theta0)
if nargin<4
    theta0 = 0;
end

if isempty(Z)
    theta = 0;
    return;
end

% Find a valid initial value for theta
max_iter = 100; iter = 0;
while ~isfinite(calcLogLikelihood(X,y,Z,theta0))
    theta0 = theta0 / 2;
    iter = iter + 1;
    assert(iter<max_iter,'Could not find valid initial value of theta for optimization');
end

warning('off','MATLAB:nearlySingularMatrix');

opts = optimoptions('fminunc','Display','off','Algorithm','Quasi-Newton');
theta = fminunc( @(x) calcLogLikelihood(X,y,Z,x) , theta0 , opts );

warning('on','MATLAB:nearlySingularMatrix');
end

%% Calculate log-likelihood associated with a given value of theta
function LL = calcLogLikelihood(X,y,Z,theta)
if nargin<4, theta = 0; end
[~,~,~,LL] = solveLME(X,y,Z,theta);
LL = nanmean(-LL);
end

%% Solve the linear mixed effects model
function [beta,bHat,covb,PLogLik,sigma2] = solveLME(X,Y,Z,theta,weights)
if nargin<5, weights = speye(size(X,1)); end
if nargin<4, theta = 0; end
if nargin<3, Z = []; end

[nT,nY] = size(Y);
nX = size(X,2);
nZ = size(Z,2);

X0 = X;
X = weights * X;
Y = weights * Y;
Z = weights * Z;

Lambda = sqrt(exp(theta)) * speye(nZ); % Isotropic covariance pattern
Iq = spdiags(ones(nZ,1),0,nZ,nZ);

[R,~,S] = chol( Lambda'*sparse(Z'*Z)*Lambda + Iq );
Q1 = ((X'*Z*Lambda)*S) / R;
R1R1t = X'*X - Q1*Q1';
R1 = cholSafe(R1R1t,'lower');

% Parameter estimates
cDeltab = R' \ (S'*((Lambda'*Z'*Y)));
cbeta = R1 \ (X'*Y - Q1*cDeltab);
beta = R1' \ cbeta;
Deltab = S*(R \ (cDeltab - Q1'*beta));
bHat = Lambda * Deltab;

% Estimate error
resid=(Y - X*beta - Z*bHat);
r2 = sum(Deltab.^2) + sum(resid.^2);

% Calculate coefficient covariance
dfe = max(nT-nX,0);
if isempty(Z)
    sigma2 = r2/dfe;
else
    sigma2 = r2/nT;
end
xr = rank(X);
invR1 = R1(1:xr,1:xr) \ eye(xr);
covb = sigma2*(invR1'*invR1);

% Calculate log-likelihood
PLogLik = (-nT/2)*( 1 + log( 2*pi*r2/nT ) ) - logDet(R);

end

function d = logDet(M)
% Safely compute the log of the determinant
% Matlab's det uses LU which is inaccurate for large matrices, and the
% determinant can easily overflow a double. Here we use QR and avoid overflowing
% by moving log to the inside, since log(a*b) = log(a)+log(b)
[~,R] = qr(M,0);
d = sum(log(diag(R)));
end

%% Cholesky decomposition that wont error if unstable
function [R,p] = cholSafe(d,varargin)

delta = eps(class(d));
I = eye(size(d));

p=1; iter=1; max_iter=1000;
while p~=0
    try
        [R,p] = chol(d + delta*I,varargin{:});
    end
    delta = 2*delta;
    iter = iter + 1;
    assert(iter<max_iter,'Could not perform cholesky factorization');
end 

end

function s = robustSigma(X,Y,Z,adj,num_params,tune,beta,bHat)
% Get robust estimate of sigma
[nT,nX] = size(X);
resid = (Y - X*beta - Z*bHat) .* adj;
[resid_s,s] = studentizeResiduals( resid , num_params );

r = resid_s/tune;
r1 = resid_s/tune - 1e-4;
r2 = resid_s/tune + 1e-4;
w = bisquare( r , 1 );
w1 = bisquare( r1 , 1 );
w2 = bisquare( r2 , 1 );
dw = (r2.*w2.^2 - r1.*w1.^2) ./ 2e-4;

a = mean(dw);
h = (1./adj).^2;
b = sum(h.*(r.*w.^2).^2)/(nT-nX);
K = 1 + (nX/nT) * (1-a) / a;
s = K*sqrt(b) * s * tune / a ;
end

function w = bisquare(resid,tune)
if isempty(tune), tune = 4.685; end
r = resid/tune;
w = (1 - r.^2) .* (r < 1 & r > -1);
end

function [resid_s,sigma] = studentizeResiduals( resid , num_params )
sa_resid = sort(abs(resid));
sigma = nanmedian(sa_resid(num_params:end)) / 0.6745;
resid_s = resid ./ sigma;
end
