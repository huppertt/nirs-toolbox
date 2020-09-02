function [b,dev,stats] = mnrfit(x,y,varargin)
%MNRFIT Fit a nominal or ordinal multinomial regression model.
%   B = MNRFIT(X,Y) fits a nominal multinomial logistic regression model for
%   the response Y and predictor matrix X.  X is an N-by-P design matrix with
%   N observations on P predictor variables.  Y is an N-by-K matrix, where
%   Y(I,J) is the number of outcomes of the multinomial category J for the
%   predictor combinations given by X(I,:).  The sample sizes for each
%   observation (rows of X and Y) are given by the row sums SUM(Y,2).
%   Alternatively, Y can be an N element column vector of scalar integers from
%   1 to K indicating the value of the response for each observation, and all
%   sample sizes are taken to be 1.  MNRFIT automatically includes intercept
%   (constant) terms; do not enter a column of ones directly into X.
%
%   The result B is a (P+1)-by-(K-1) matrix of estimates, where each column
%   corresponds to the estimated intercept term and predictor coefficients,
%   one for each of the first (K-1) multinomial categories.  The estimates for
%   the K-th category are taken to be zero.
%
%   MNRFIT treats NaNs in X and Y as missing data, and removes the
%   corresponding observations.
%
%   B = MNRFIT(X,Y,'PARAM1',val1,'PARAM2',val2,...) allows you to
%   specify optional parameter name/value pairs to control the model fit.
%   Parameters are:
%
%      'model' - the type of model to fit, one of the text strings 'nominal'
%         (the default), 'ordinal', or 'hierarchical'.
%
%      'interactions' - determines whether the model includes an interaction
%         between the multinomial categories and the coefficients.  Specify as
%         'off' to fit a model with a common set of coefficients for the
%         predictor variables, across all multinomial categories.  This is
%         often described as "parallel regression".  Specify as 'on' to fit a
%         model with different coefficients across categories.  In all cases,
%         the model has different intercepts across categories.  Thus, B is a
%         vector containing K-1+P coefficient estimates when 'interaction' is
%         'off', and a (P+1)-by-(K-1) matrix when it is 'on'. The default is
%         'off' for ordinal models, and 'on' for nominal and hierarchical
%         models.
%
%      'link' - the link function to use for ordinal and hierarchical models.
%         The link function defines the relationship g(mu_ij) = x_i*b_j
%         between the mean response for the i-th observation in the j-th
%         category, mu_ij, and the linear combination of predictors x_i*b_j.
%         Specify the link parameter value as one of the text strings 'logit'
%         (the default), 'probit', 'comploglog', or 'loglog'.  You may not
%         specify the 'link' parameter for nominal models; these always use a
%         multivariate logistic link.
%
%      'estdisp' - specify as 'on' to estimate a dispersion parameter for
%         the multinomial distribution in computing standard errors, or 'off'
%         (the default) to use the theoretical dispersion value of 1.
%
%   [B,DEV] = MNRFIT(...) returns the deviance of the fit.
%
%   [B,DEV,STATS] = MNRFIT(...) returns a structure that contains the
%   following fields:
%       'dfe'       degrees of freedom for error
%       's'         theoretical or estimated dispersion parameter
%       'sfit'      estimated dispersion parameter
%       'se'        standard errors of coefficient estimates B
%       'coeffcorr' correlation matrix for B
%       'covb'      estimated covariance matrix for B
%       't'         t statistics for B
%       'p'         p-values for B
%       'resid'     residuals
%       'residp'    Pearson residuals
%       'residd'    deviance residuals
%
%   Example:  Fit multinomial logistic regression models to data with one
%   predictor variable and three categories in the response variable.
%
%     x = [-3 -2 -1 0 1 2 3]';
%     Y = [1 11 13; 2 9 14; 6 14 5; 5 10 10; 5 14 6; 7 13 5; 8 11 6];
%     bar(x,Y,'stacked'); ylim([0 25]);
%
%     % Fit a nominal model for the individual response category probabilities,
%     % with separate slopes on the single predictor variable, x, for each
%     % category.  The first row of betaHatNom contains the intercept terms for
%     % the first two response categories.  The second row contains the slopes.
%     betaHatNom = mnrfit(x,Y,'model','nominal','interactions','on')
%
%     % Compute the predicted probabilities for the three response categories.
%     xx = linspace(-4,4)';
%     pHatNom = mnrval(betaHatNom,xx,'model','nominal','interactions','on');
%     line(xx,cumsum(25*pHatNom,2),'LineWidth',2);
%
%     % Fit a "parallel" ordinal model for the cumulative response category
%     % probabilities, with a common slope on the single predictor variable, x,
%     % across all categories.  The first two elements of betaHatOrd are the
%     % intercept terms for the first two response categories.  The last element
%     % of betaHatOrd is the common slope.
%     betaHatOrd = mnrfit(x,Y,'model','ordinal','interactions','off')
%
%     % Compute the predicted cumulative probabilities for the first two response
%     % categories.  The cumulative probability for the third category is always 1.
%     pHatOrd = mnrval(betaHatOrd,xx,'type','cumulative','model','ordinal','interactions','off');
%     bar(x,cumsum(Y,2),'grouped'); ylim([0 25]);
%     line(xx,25*pHatOrd,'LineWidth',2);
%
%   See also MNRVAL, GLMFIT, GLMVAL, REGRESS, REGSTATS.

%   References:
%      [1] McCullagh, P., and J.A. Nelder (1990) Generalized Linear
%          Models, 2nd edition, Chapman&Hall/CRC Press.

%   Copyright 2006-2018 The MathWorks, Inc.

dorobust=true;  % TJH- flag for robust fitting


if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

narginchk(2,Inf);

pnames = {  'model' 'interactions' 'link' 'estdisp'};
dflts =  {'nominal'             []     []     'off'};
[model,interactions,link,estdisp] = ...
    internal.stats.parseArgs(pnames, dflts, varargin{:});
iterlim = 100;
tolpos = eps(class(x))^(3/4);
model = internal.stats.getParamVal(model,{'nominal','ordinal','hierarchical'},'''model''');

if isempty(interactions)
    % Default is 'off' for ordinal models, 'on' for nominal or hierarchical
    parallel = strcmp(model,'ordinal');
elseif isequal(interactions,'on')
    parallel = false;
elseif isequal(interactions,'off')
    parallel = true;
elseif islogical(interactions)
    parallel = ~interactions;
else % ~islogical(interactions)
    error(message('stats:mnrfit:BadInteractions'));
end
if parallel && strcmp(model,'nominal')
    % A nominal model with no interactions is the same as having no predictors.
    warning(message('stats:mnrfit:NominalNoInteractions'));
    x = zeros(size(x,1),0,'like',x);
end

% Categorical responses
if isa(y,'categorical')
    y = grp2idx(y);
end

dataClass = superiorfloat(x,y);

if isempty(link)
    link = 'logit';
elseif ~isempty(link) && strcmp(model,'nominal')
    error(message('stats:mnrfit:LinkNotAllowed'));
else
    link = internal.stats.getParamVal(link,{'logit' 'probit' 'comploglog' 'loglog'},'''link''');
end
[flink,dlink,ilink] = stattestlink(link,dataClass);

if isequal(estdisp,'on')
    estdisp = true;
elseif isequal(estdisp,'off')
    estdisp = false;
elseif ~islogical(estdisp)
    error(message('stats:mnrfit:BadEstDisp'));
end

% Remove missing values from the data.  Also turns row vectors into columns.
[anybad,wasnan,y,x] = statremovenan(y,x);
if anybad
    error(message('stats:mnrfit:InputSizeMismatch'));
end
p = size(x,2);
[n,k] = size(y);
if n == 0
    error(message('stats:mnrfit:NoData'));
end

if k == 1
    if min(y) < 1 || any(y ~= floor(y))
        error(message('stats:mnrfit:BadY'));
    end
    y = accumarray({(1:n)' y},ones(dataClass));
    k = size(y,2);
    m = ones(n,1,dataClass);
else
    m = sum(y,2);
end
if parallel
    pstar = k - 1 + p;
    dfe = n * (k-1) - pstar;
else
    pstar = p + 1;
    dfe = (n-pstar) * (k-1);
end
    if strcmp(model,'hierarchical')
        if nargout < 3
            [b,dev] = hierarchicalFit(x,y,m,link,n,k,p,pstar,parallel,estdisp);
        else
            [b,dev,stats] = ...
                hierarchicalFit(x,y,m,link,n,k,p,pstar,parallel,estdisp);
        end
    else
        % Set up initial estimates from the data themselves
        pi = (y+0.5) ./ (m+k/2); % shrink raw proportions toward uniform
        if strcmp(model,'nominal')
            [b,hess,pi] = nominalFit(x,y,m,pi,n,k,p,pstar,parallel,iterlim,tolpos);
        else % 'ordinal'
            z = cumsum(y(:,1:(k-1)),2);
            [b,hess,pi,gam] = ...
                ordinalFit(x,z,m,pi,flink,ilink,dlink,n,k,p,pstar,parallel,iterlim,tolpos);
        end
        
        % Deviance residuals - one for each vector observation of cell counts
        mu = pi .* m;
        D = zeros(size(y),dataClass);
        t = (y > 0); % avoid 0*log(0), but let (pi==0) & (y>0) happen
        D(t) = 2 * y(t) .* log(y(t) ./ mu(t));
        rd = sum(D,2);
        dev = sum(rd);
        
        if nargout > 2
            % The Pearson residuals in terms of y and pi are not equivalent to
            % those computed using z and gamma.  Use the appropriate version to
            % estimate dispersion.
            if strcmp(model,'nominal')
                r = y - pi .* m;
                pi = constrain(pi,tolpos,1-tolpos); % don't allow Inf variance
                rp = r ./ sqrt(pi .* (1 - pi) .* m);
                sigsq = ((k-1)/k) * sum(sum(rp .* rp)) ./ dfe; % bias corrected
            elseif strcmp(model,'ordinal')
                r = z - gam .* m;
                rp = r ./ sqrt(gam .* (1 - gam) .* m);
                sigsq = sum(sum(rp .* rp)) ./ dfe;
            end
            stats.beta = b;
            stats.dfe = dfe;
            if dfe > 0
                stats.sfit = sqrt(sigsq);
            else
                stats.sfit = NaN;
            end
            if estdisp
                stats.s = stats.sfit;
                rp = rp ./ stats.sfit;
            else
                stats.s = ones(dataClass);
            end
            stats.estdisp = estdisp;
            
            if ~isnan(stats.s) % dfe > 0 or estdisp == 'off'
                % bcov = inv(hess); bcov = (bcov + bcov')/2;
                bcov = linsolve(hess,eye(size(hess)),struct('SYM',true,'POSDEF',true));
                if estdisp
                    bcov = bcov * sigsq;
                end
                se = sqrt(diag(bcov));
                stats.covb = bcov;
                stats.coeffcorr = bcov ./ (se*se');
                if ~parallel
                    se = reshape(se,pstar,k-1);
                end
                stats.se = se;
                stats.t = b ./ se;
                if estdisp
                    stats.p = 2 * tcdf(-abs(stats.t), dfe);
                else
                    stats.p = 2 * normcdf(-abs(stats.t));
                end
            else
                stats.se = NaN(size(b),dataClass);
                stats.coeffcorr = NaN(numel(b),dataClass);
                stats.t = NaN(size(b),dataClass);
                stats.p = NaN(size(b),dataClass);
            end
            stats.resid = r;
            stats.residp = rp;
            stats.residd = rd;
        end
    end
   

 if nargout > 2 && any(wasnan)
        stats.resid  = statinsertnan(wasnan, stats.resid);
        stats.residp = statinsertnan(wasnan, stats.residp);
        stats.residd = statinsertnan(wasnan,stats.residd);
    end

%------------------------------------------------------------------------
function [b,XWX,pi,gam] = ordinalFit(x,z,m,pi,flink,ilink,dlink,n,k,p,pstar,parallel,iterLim,tolpos)

kron1 = repmat(1:k-1,pstar,1);
kron2 = repmat((1:pstar)',1,k-1);

gam = cumsum(pi(:,1:(k-1)),2);
eta = flink(gam);

% Main IRLS loop
iter = 0;
seps = sqrt(eps); % don't depend on class
convcrit = 1e-6;
b = 0;
while iter <= iterLim
    iter = iter + 1;
    
    % d.gamma(i,)/d.eta(i,) is actually (k-1) by (k-1) but diagonal,
    % so can store d.mu/d.eta as n by (k-1) even though it is really
    % n by (k-1) by (k-1)
    mu = m .* gam;
    deta = dlink(gam) ./ m; % d(eta)/d(mu)
    dmu = 1 ./ deta;  % d(mu)/d(eta)
    
    % Adjusted dependent variate
    Z = eta + deta.*(z - mu);
    
    % Tridiagonal symmetric weight matrix (scaled by m)
    diagW = dmu .* dmu .* (1./pi(:,1:(k-1)) + 1./pi(:,2:k));
    offdiagW = -(dmu(:,1:(k-2)) .* dmu(:,2:k-1)) ./ pi(:,2:(k-1));
    
    % Update the coefficient estimates.
    b_old = b;
    XWX = 0;
    XWZ = 0;
    for i = 1:n
        W = (1./m(i)) .* (diag(diagW(i,:)) + ...
            diag(offdiagW(i,:),1) + diag(offdiagW(i,:),-1));
        if p > 0
            % The first step for a nonparallel model can be wild, so fit
            % a parallel model for the first iteration, regardless
            if parallel || (iter==1)
                % Do these computations, but more efficiently
                % Xstar = [eye(k-1) repmat(x(i,:),k-1,1)];
                % XWX = XWX + Xstar'*W*Xstar;
                % XWZ = XWZ + Xstar'*W*Z(i,:)';
                xi = x(i,:);
                OneW = sum(W,1);
                xOneW = xi'*OneW;
                XWX = XWX + [W      xOneW'; ...
                    xOneW  sum(OneW)*(xi'*xi)];
                XWZ = XWZ + [W; xOneW] * Z(i,:)';
            else
                xstar = [1 x(i,:)];
                % Do these computations, but more efficiently
                % XWX = XWX + kron(W, xstar'*xstar);
                % XWZ = XWZ + kron(W*Z(i,:)', xstar');
                XWX = XWX + W(kron1,kron1) .* (xstar(1,kron2)'*xstar(1,kron2));
                WZ = Z(i,:)*W;
                XWZ = XWZ + WZ(1,kron1)' .* xstar(1,kron2)';
            end
        else
            XWX = XWX + W;
            XWZ = XWZ + W * Z(i,:)';
        end
    end
    b = XWX \ XWZ;
    
    % Update the linear predictors.
    eta_old = eta;
    if parallel
        if p > 0
            eta = repmat(b(1:(k-1))',n,1) + repmat(x*b(k:pstar),1,k-1);
        else
            eta = repmat(b',n,1);
        end
    else
        if iter == 1
            % the first iteration was a parallel fit, transform those
            % estimates to the equivalent non-parallel format.
            b = [b(1:k-1)'; repmat(b(k:end),1,k-1)];
        else
            % Convert from vector to the matrix format.
            b = reshape(b,pstar,k-1);
        end
        if p > 0
            eta = b(1,:) + x*b(2:pstar,:); % row plus full matrix
        else
            eta = repmat(b,n,1);
        end
    end
    
    gam = ilink(eta);
    diffgam = diff(gam,[],2);
    pi = [gam(:,1) diffgam 1-gam(:,k-1)];
    
    % Check stopping conditions
    cvgTest = abs(b-b_old) > convcrit * max(seps, abs(b_old));
    if iter<=iterLim && any(cvgTest(:))
        % Except at final iteration, try to avoid a bad step
        for backstep = 0:10
            if ~any(pi(:)<0)
                break
            end
            
            % Try a shorter step in the same direction if we have negative
            % probabilities. eta_old is feasible, even on iteration 1,
            % though it may be that no coefficients will give it.
            if backstep < 10
                eta = eta_old + (eta - eta_old)/5;
                gam = ilink(eta);
                diffgam = diff(gam,[],2);
                pi = [gam(:,1) diffgam 1-gam(:,k-1)];
            end
        end
    end
    
    % If some probabilities are zero, that could cause problems in fitting,
    % so force them to a small positive number and update gam to match.
    if ~all(pi(:) > tolpos)
        pi = max(pi,tolpos);
        pi = pi ./ sum(pi,2); % force rows to sum to 1
        gam = cumsum(pi(:,1:k-1),2);
        eta = flink(gam);
    end
    
    % Check stopping conditions
    if (~any(cvgTest(:))), break; end
end
if iter > iterLim
    warning(message('stats:mnrfit:IterOrEvalLimit'));
end


%------------------------------------------------------------------------
function [b,XWX,pi] = nominalFit(x,y,m,pi,n,k,p,pstar,parallel,iterLim,tolpos)

kron1 = repmat(1:k-1,pstar,1);
kron2 = repmat((1:pstar)',1,k-1);

eta = log(pi);

% Main IRLS loop
iter = 0;
seps = sqrt(eps); % don't depend on class
convcrit = 1e-6;
b = 0;
dataClass = superiorfloat(x,y);
lowerBnd = log(eps(dataClass));
upperBnd = -lowerBnd;
refine = false; % flag to indicate solution may need more refinement


while iter <= iterLim
    iter = iter + 1;
    mu = m .* pi;
    
    % Updated the coefficient estimates.
    b_old = b;
    XWX = 0;
    XWZ = 0;
    
    
    for i = 1:n
        W = diag(mu(i,:)) - mu(i,:)'*pi(i,:);
        
        % Adjusted dependent variate
        Z = eta(i,:)*W + (y(i,:) - mu(i,:));
        
        if p > 0 % parallel models with p>0 have been weeded out
            xstar = [1 x(i,:)];
            % Do these computations, but more efficiently
            % XWX = XWX + kron(W(1:k-1,1:k-1), xstar'*xstar);
            % XWZ = XWZ + kron(Z(1:k-1)', xstar');
            XWX = XWX + W(kron1,kron1) .* (xstar(1,kron2)'*xstar(1,kron2));
            XWZ = XWZ + Z(1,kron1)' .* xstar(1,kron2)';
        else
            XWX = XWX + W(1:k-1,1:k-1);
            XWZ = XWZ + Z(1:k-1)';
        end
    end
    
    b = XWX \ XWZ;
    
    
    % Update the linear predictors.
    eta_old = eta;
    if parallel % parallel models with p>0 have been simplified already
        eta = repmat(b',n,1);
    else
        b = reshape(b,pstar,k-1);
        if p > 0
            eta = b(1,:) + x*b(2:pstar,:); % row plus full matrix
        else
            eta = repmat(b,n,1);
        end
    end
    eta = [eta zeros(n,1,'like',eta)];
    
    % Check stopping conditions
    cvgTest = abs(b-b_old) > convcrit * max(seps, abs(b_old));
    if iter>iterLim || ~any(cvgTest(:))
        pi = exp(eta-max(eta,[],2));   % rescale but do not constrain
        pi = pi ./ sum(pi,2);          % make rows sum to 1
        break
    end
    
    refine = false;
    for backstep = 0:10
        % Update the predicted category probabilities; constrain to a
        % reasonable range to avoid problems during fitting.
        pi = exp(constrain(eta-max(eta,[],2),lowerBnd,upperBnd));
        pi = pi ./ sum(pi,2); % make rows sum to 1
        
        % If all observations have positive category probabilities,
        % we can take the step as is.
        if all(pi(:) > tolpos)
            break;
            
            % In either of the following cases the coefficients vector will not
            % yield the pi and eta values, but we have not converged so we will
            % make them consistent in the next iteration. We will not test for
            % convergence while backtracking, because we don't have coefficient
            % estimates that give the shorter steps.
            
            % Otherwise try a shorter step in the same direction.  eta_old is
            % feasible, even on the first iteration.
        elseif backstep < 10
            eta = eta_old + (eta - eta_old)/5;
            refine = true;
            
            % If the step direction just isn't working out, force the
            % category probabilities to be positive, and make the linear
            % predictors compatible with that.
        else
            pi = max(pi,tolpos);
            pi = pi ./ sum(pi,2); % make rows sum to 1
            eta = log(pi);
            refine = true;
            break;
        end
    end
end
if refine && ~parallel && p>0
    % The iterative reweighted least squares procedure may stall as fitted
    % probabilities get close to 0 or 1. Try to refine the results with
    % some iterations of fminsearch. Must not have single data type.
    [b,pi] = refineEstimate(b,x,y,pi,pstar);
end
if iter > iterLim
    warning(message('stats:mnrfit:IterOrEvalLimit'));
end

function [b,pi] = refineEstimate(b,x,y,pi,pstar)
% try to refine estimate; if any failure just leave as is
try
    op = optimset('Display','none');
    newb = fminsearch(@(newb)calcdev(newb,size(b),x,y),double(b),op);
    b = reshape(newb,size(b));
    eta = b(1,:) + x*b(2:pstar,:); % row plus full matrix
    n = size(eta,1);
    eta = [eta zeros(n,1,'like',eta)];
    
    pi = exp(eta-max(eta,[],2));   % rescale but do not constrain
    pi = pi ./ sum(pi,2);          % make rows sum to 1
catch
end

function dev = calcdev(B,sz,x,y) % deviance calculation for fminsearch
B = reshape(B,sz);
p = mnrval(B,x);
t = y>0;
dev = -2*sum(y(t).*log(p(t)));

%------------------------------------------------------------------------
function [b,dev,stats] = hierarchicalFit(x,y,m,link,n,k,~,pstar,parallel,estdisp)

dataClass = superiorfloat(x,y);

% Compute the sample sizes for the conditional binomial observations.  Some
% might be zero, rely on glmfit to ignore those, tell us the right dfe, and
% return NaN residuals there.
m = [m repmat(m,1,k-2)-cumsum(y(:,1:(k-2)),2)];

warnStateSaved = warning('off','stats:glmfit:IterationLimit');
[wmsgSaved,widSaved] = lastwarn;
lastwarn(''); % clear this so we can look for a new iter limit warning
needToWarn = false;
try
    if parallel
        % Same slopes for the categories, fit a single binomial model by
        % transforming the multinomial observations into conditional binomial
        % observations.
        ii = repmat(1:n,1,k-1);
        jj = repmat(1:k-1,n,1);
        dummyvars = eye(k-1,k-1,dataClass);
        xstar = [dummyvars(jj,:) x(ii,:)];
        ystar = y(:,1:k-1);
        if estdisp
            estdisp = 'on';
        else
            estdisp = 'off';
        end
        if nargout < 3
            [b,dev] = glmfit(xstar,[ystar(:) m(:)],'binomial',...
                'link',link,'constant','off','estdisp',estdisp);
            needToWarn = checkForIterWarn(needToWarn);
        else
            [b,dev,stats] = glmfit(xstar,[ystar(:) m(:)],'binomial', ...
                'link',link,'constant','off','estdisp',estdisp);
            needToWarn = checkForIterWarn(needToWarn);
            stats.resid = reshape(stats.resid,n,k-1);
            stats.residp = reshape(stats.residp,n,k-1);
            stats.residd = sum(reshape(stats.residd,n,k-1),2);
            stats = rmfield(stats,'resida');
        end
        
    else % ~parallel
        % Separate slopes for the categories, fit a sequence of conditional
        % binomial models
        b = zeros(pstar,k-1,dataClass);
        dev = zeros(dataClass);
        if nargout < 3
            for j = 1:k-1
                [b(:,j),d] = glmfit(x,[y(:,j) m(:,j)], 'binomial','link',link);
                needToWarn = checkForIterWarn(needToWarn);
                dev = dev + d;
            end
        else
            stats = struct('beta',zeros(pstar,k-1,dataClass), ...
                'dfe',zeros(dataClass), ...
                'sfit',NaN(dataClass), ...
                's',ones(dataClass), ...
                'estdisp',estdisp, ...
                'se',zeros(pstar,k-1,dataClass), ...
                'coeffcorr',zeros(pstar*(k-1),dataClass), ...
                't',zeros(pstar,k-1,dataClass), ...
                'p',zeros(pstar,k-1,dataClass), ...
                'resid',zeros(n,k-1,dataClass), ...
                'residp',zeros(n,k-1,dataClass), ...
                'residd',zeros(n,1,dataClass));
            for j = 1:k-1
                [b(:,j),d,s] = glmfit(x,[y(:,j) m(:,j)], 'binomial','link',link);
                needToWarn = checkForIterWarn(needToWarn);
                dev = dev + d;
                stats.beta(:,j) = b(:,j);
                stats.dfe = stats.dfe + s.dfe; % not n-pstar if some m's are zero
                stats.se(:,j) = s.se;
                jj = (j-1)*pstar + (1:pstar);
                stats.coeffcorr(jj,jj) = s.coeffcorr;
                stats.p(:,j) = s.p;
                stats.t(:,j) = s.t;
                stats.resid(:,j)  = s.resid;
                stats.residp(:,j) = s.residp;
                stats.residd = stats.residd + s.residd;
            end
            if stats.dfe > 0
                % Weed out the NaN residuals caused by zero conditional sizes
                % when computing dispersion.
                t = ~isnan(stats.residp(:));
                sigsq = sum(stats.residp(t) .* stats.residp(t)) ./ stats.dfe;
                stats.sfit = sqrt(sigsq);
            else
                % stats.sfit already NaN
            end
            if estdisp
                sigma = stats.sfit;
                stats.s = sigma;
                stats.residp = stats.residp ./ sigma;
                stats.se = stats.se .* sigma;
                stats.t = stats.t ./ sigma;
                stats.p = 2 * tcdf(-abs(stats.t), stats.dfe);
            else
                % stats.s already 1
            end
        end
    end
catch ME
    warning(warnStateSaved);
    rethrow(ME);
end
[~,wid] = lastwarn;
if needToWarn
    warning(message('stats:mnrfit:IterOrEvalLimit'));
elseif ~isempty(widSaved) && isempty(wid)
    % Restore any pre-existing warning if there was not a new one.
    lastwarn(wmsgSaved,widSaved);
end
warning(warnStateSaved);

function needToWarn = checkForIterWarn(needToWarn)
[~,wid] = lastwarn;
needToWarn = needToWarn || strcmp(wid,'stats:glmfit:IterationLimit');

function x = constrain(x,lower,upper)
% Constrain between upper and lower limits. Avoid max(x,lower) because we
% want to retain NaN values when they occur.
x(x<lower) = lower;
x(x>upper) = upper;
