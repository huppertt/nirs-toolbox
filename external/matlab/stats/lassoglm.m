function [B,stats] = lassoglm(x,y,distr,varargin)
%LASSOGLM Perform lasso or elastic net regularization for a generalized linear model.
%   [B,STATS] = LASSOGLM(X,Y,DISTR,...) Performs L1-penalized maximum likelihood 
%   fits (lasso) relating the predictors in X to the responses in Y, or 
%   fits subject to joint L1- and L2 constraints (elastic net).
%   The default is a lasso-style fit, that is, a maximum likelihood fit
%   subject to a constraint on the L1-norm of the coefficients B.
%
%   LASSOGLM accepts all the command line parameters of the LASSO function 
%   and it accepts command line parameters of the GLMFIT function, with 
%   the following exceptions. LASSOGLM also accepts the arguments 'Link'
%   and 'Offset' of GLMFIT (they can also be lower case 'link' or 'offset').  
%   LASSOGLM It does not accept the argument 'estdisp' of GLMFIT, or the 
%   argument 'constant'.  LASSOGLM does not calculate standard errors
%   or covariances among the coefficients, as GLMFIT does.
%   LASSOGLM always calculates an Intercept term. 
%   
%   Positional parameters:
%
%     X                A numeric matrix (dimension, say, NxP)
%     Y                A numeric vector of length N.  If DISTR is
%                      'binomial', Y may be a binary vector indicating
%                      success/failure, and the total number of trials is
%                      taken to be 1 for all observations.  If DISTR is
%                      'binomial', Y may also be a two column matrix, the 
%                      first column containing the number of successes for
%                      each observation, and the second column containing
%                      the total number of trials.
%     DISTR            The distributional family for the non-systematic
%                      variation in the responses.  Acceptable values for DISTR are
%                      'normal', 'binomial', 'poisson', 'gamma', and 'inverse gaussian'.  
%                      By default, the distribution is fit using the canonical link 
%                      corresponding to DISTR.  Other link functions may be
%                      specified using the optional parameter 'link'.
%   
%   Optional input parameters:  
%
%     'Weights'        Observation weights.  Must be a vector of non-negative
%                      values, of the same length as columns of X.  At least
%                      two values must be positive. (default ones(N,1) or 
%                      equivalently (1/N)*ones(N,1)).
%     'Alpha'          Elastic net mixing value, or the relative balance
%                      between L2 and L1 penalty (default 1, range (0,1]).
%                      Alpha=1 ==> lasso, otherwise elastic net.
%                      Alpha near zero ==> nearly ridge regression.
%     'NumLambda'      The number of lambda values to use, if the parameter
%                      'Lambda' is not supplied (default 100). Ignored
%                      if 'Lambda' is supplied. LASSOGLM may return fewer
%                      fits than specified by 'NumLambda' if the deviance
%                      of the fits drops below a threshold percentage of
%                      the null deviance (deviance of the fit without
%                      any predictors X).
%     'LambdaRatio'    Ratio between the minimum value and maximum value of
%                      lambda to generate, if the  parameter "Lambda" is not 
%                      supplied.  Legal range is [0,1). Default is 0.0001.
%                      If 'LambdaRatio' is zero, LASSOGLM will generate its
%                      default sequence of lambda values but replace the
%                      smallest value in this sequence with the value zero.
%                      'LambdaRatio' is ignored if 'Lambda' is supplied.
%     'Lambda'         Lambda values. Will be returned in return argument
%                      STATS in ascending order. The default is to have
%                      LASSOGLM generate a sequence of lambda values, based 
%                      on 'NumLambda' and 'LambdaRatio'. LASSOGLM will generate 
%                      a sequence, based on the values in X and Y, such that 
%                      the largest Lambda value is estimated to be just 
%                      sufficient to produce all zero coefficients B. 
%                      You may supply a vector of real, non-negative values 
%                      of lambda for LASSOGLM to use, in place of its default
%                      sequence.  If you supply a value for 'Lambda', 
%                      'NumLambda' and 'LambdaRatio' are ignored.
%     'DFmax'          Maximum number of non-zero coefficients in the model.
%                      Can be useful with large numbers of predictors.
%                      Results only for lambda values that satisfy this
%                      degree of sparseness will be returned.  Default is
%                      to not limit the number of non-zero coefficients.
%     'Standardize'    Whether to scale X prior to fitting the model
%                      sequence. This affects whether the regularization is
%                      applied to the coefficients on the standardized
%                      scale or the original scale. The results are always
%                      presented on the original data scale. Default is
%                      TRUE, do scale X.
%     'RelTol'         Convergence threshold for coordinate descent algorithm.
%                      The coordinate descent iterations will terminate
%                      when the relative change in the size of the
%                      estimated coefficients B drops below this threshold.
%                      Default: 1e-4. Legal range is (0,1).
%     'CV'             If present, indicates the method used to compute Deviance.
%                      When 'CV' is a positive integer K, LASSOGLM uses K-fold
%                      cross-validation.  Set 'CV' to a cross-validation 
%                      partition, created using CVPARTITION, to use other
%                      forms of cross-validation. You cannot use a
%                      'Leaveout' partition with LASSOGLM.                
%                      When 'CV' is 'resubstitution', LASSOGLM uses X and Y 
%                      both to fit the model and to estimate the deviance 
%                      of the fitted model, without cross-validation.  
%                      The default is 'resubstitution'.
%     'MCReps'         A positive integer indicating the number of Monte-Carlo
%                      repetitions for cross-validation.  The default value is 1.
%                      If 'CV' is 'resubstitution' or a cvpartition of type
%                      'resubstitution', 'MCReps' must be 1.  If 'CV' is a
%                      cvpartition of type 'holdout', then 'MCReps' must be
%                      greater than one.
%     'PredictorNames' A cell array of names for the predictor variables,
%                      in the order in which they appear in X. 
%                      Default: {}
%     'Options'        A structure that contains options specifying whether to
%                      conduct cross-validation evaluations in parallel, and
%                      options specifying how to use random numbers when computing
%                      cross validation partitions. This argument can be created
%                      by a call to STATSET. CROSSVAL uses the following fields:
%                        'UseParallel'
%                        'UseSubstreams'
%                        'Streams'
%                      For information on these fields see PARALLELSTATS.
%                      NOTE: If supplied, 'Streams' must be of length one.
%     'Link'           The link function to use in place of the canonical link.
%                      The link function defines the relationship f(mu) = x*b
%                      between the mean response mu and the linear combination of
%                      predictors x*b.  Specify the link parameter value as one of
%                      - the text strings 'identity', 'log', 'logit', 'probit',
%                        'comploglog', 'reciprocal', 'loglog', or
%                      - an exponent P defining the power link, mu = (x*b)^P 
%                        for x*b > 0, or
%                      - a cell array of the form {FL FD FI}, containing three
%                        function handles, created using @, that define the link (FL),
%                        the derivative of the link (FD), and the inverse link (FI).  
%     'Offset'         A vector to use as an additional predictor variable, but
%                      with a coefficient value fixed at 1.0.  Default is
%                      not to utilize an offset variable.
% 
%   Return values:
%     B                The fitted coefficients for each model. 
%                      B will have dimension PxL, where 
%                      P = size(X,2) is the number of predictors, and
%                      L = length(STATS.Lambda).
%     STATS            STATS is a struct that contains information about the
%                      sequence of model fits corresponding to the columns
%                      of B. STATS contains the following fields:
%
%       'Intercept'    The intercept term for each model. Dimension 1xL.
%       'Lambda'       The sequence of lambda penalties used, in ascending order. 
%                      Dimension 1xL.
%       'Alpha'        The elastic net mixing value that was used.
%       'DF'           The number of nonzero coefficients in B for each
%                      value of lambda. Dimension 1xL.
%       'Deviance'     The deviance of the fitted model for each value of
%                      lambda. If cross-validation was performed, the values 
%                      for 'Deviance' represent the estimated expected 
%                      deviance of the model applied to new data, as 
%                      calculated by cross-validation. Otherwise, 
%                      'Deviance' is the deviance of the fitted model 
%                      applied to the data used to perform the fit. 
%                      Dimension 1xL.
%
%     If cross-validation was performed, STATS also includes the following
%     fields:
%
%       'SE'                The standard error of 'Deviance' for each lambda, as
%                           calculated during cross-validation. Dimension 1xL.
%       'LambdaMinDeviance' The lambda value with minimum expected deviance, as 
%                           calculated during cross-validation. Scalar.
%       'Lambda1SE'         The largest lambda such that 'Deviance' is within 
%                           one standard error of the minimum. Scalar.
%       'IndexMinDeviance'  The index of Lambda with value LambdaMinMSE. Scalar.
%       'Index1SE'          The index of Lambda with value Lambda1SE. Scalar.
%
%   See also lassoPlot, lasso, ridge, parallelstats, glmfit.

%   References: 
%   [1] Tibshirani, R. (1996) Regression shrinkage and selection
%       via the lasso. Journal of the Royal Statistical Society,
%       Series B, Vol 58, No. 1, pp. 267-288.
%   [2] Zou, H. and T. Hastie. (2005) Regularization and variable
%       selection via the elastic net. Journal of the Royal Statistical
%       Society, Series B, Vol. 67, No. 2, pp. 301-320.
%   [3] Friedman, J., R. Tibshirani, and T. Hastie. (2010) Regularization
%       paths for generalized linear models via coordinate descent.
%       Journal of Statistical Software, Vol 33, No. 1,
%       http://www.jstatsoft.org/v33/i01.
%   [4] Hastie, T., R. Tibshirani, and J. Friedman. (2008) The Elements
%       of Statistical Learning, 2nd edition, Springer, New York.
%   [5] Dobson, A.J. (2002) An Introduction to Generalized Linear
%       Models, 2nd edition, Chapman&Hall/CRC Press.
%   [6] McCullagh, P., and J.A. Nelder (1989) Generalized Linear
%       Models, 2nd edition, Chapman&Hall/CRC Press.
%   [7] Collett, D. (2003) Modelling Binary Data, 2nd edition,
%       Chapman&Hall/CRC Press.

%   Copyright 2011-2014 The MathWorks, Inc.


if nargin < 2
    error(message('stats:lassoGlm:TooFewInputs'));
end

if nargin < 3 || isempty(distr), distr = 'normal'; end

paramNames = {     'link' 'offset' 'weights'};
paramDflts = {'canonical'  []       []};
[link,offset,pwts,~,varargin] = ...
                    internal.stats.parseArgs(paramNames, paramDflts, varargin{:});
                
% Read in the optional parameters pertinent to regularization (eg, lasso)
LRdefault = 1e-4;
pnames = { 'alpha' 'numlambda' 'lambdaratio' 'lambda' ...
    'dfmax' 'standardize' 'reltol' 'cv' 'mcreps' ...
    'predictornames' 'options' };
dflts  = {  1       100       LRdefault     []      ...
     []      true          1e-4    'resubstitution'  1 ...
     {}               []};
[alpha, nLambda, lambdaRatio, lambda, ...
    dfmax, standardize, reltol, cvp, mcreps, predictorNames, parallelOptions] ...
     = internal.stats.parseArgs(pnames, dflts, varargin{:});
 
if ~isempty(lambda)
    userSuppliedLambda = true;
else
    userSuppliedLambda = false;
end

% X a real 2D matrix
if ~ismatrix(x) || length(size(x)) ~= 2 || ~isreal(x)
    error(message('stats:lassoGlm:XnotaReal2DMatrix'));
end

% We need at least two observations.
if isempty(x) || size(x,1) < 2
    error(message('stats:lassoGlm:TooFewObservations'));
end

% Categorical responses 'binomial'
if isa(y,'categorical')
    [y, classname] = grp2idx(y); 
    nc = length(classname);
    if nc > 2
        error(message('stats:glmfit:TwoLevelCategory'));
    end
    y(y==1) = 0;
    y(y==2) = 1;
end

% Number of Predictors
P = size(x,2);

% Head off potential cruft in the command window.
wsIllConditioned2 = warning('off','stats:glmfit:IllConditioned');
cleanupIllConditioned2 = onCleanup(@() warning(wsIllConditioned2));

% Sanity checking on predictors, responses, weights and offset parameter,
% and removal of NaNs and Infs from same. Also, conversion of the
% two-column form of binomial response to a proportion.
[X, Y, offset, pwts, dataClass, nTrials, binomialTwoColumn] = ...
    glmProcessData(x, y, distr, 'off', offset, pwts);

[~,sqrtvarFun,devFun,linkFun,dlinkFun,ilinkFun,link,mu,eta,muLims,isCanonical,dlinkFunCanonical] = ...
    glmProcessDistrAndLink(Y,distr,link,'off',nTrials,dataClass);

[X,Y,pwts,nLambda,lambda,dfmax,cvp,mcreps,predictorNames,ever_active] = ...
    processLassoParameters(X,Y,pwts,alpha,nLambda,lambdaRatio,lambda,dfmax, ...
    standardize,reltol,cvp,mcreps,predictorNames);

% Compute the amount of penalty at which all coefficients shrink to zero.
[lambdaMax, nullDev, nullIntercept] = computeLambdaMax(X, Y, pwts, alpha, standardize, ...
    distr, link, dlinkFun, offset, isCanonical, dlinkFunCanonical, devFun);

% If the command line did not provide a sequence of penalty terms,
% generate a sequence.
if isempty(lambda)
    lambda = computeLambdaSequence(lambdaMax, nLambda, lambdaRatio, LRdefault);
end

% Also, to help control saturated fits with binomial outcomes,
% also alter the glmfit values for muLims
if strcmp(distr,'binomial')
    muLims = [1.0e-5 1.0-1.0e-5];
end

% Will convert from lambda penalty matched to average log-likelihood
% per observation (which is what lasso uses) to an equivalent
% penalty for total log-likelihood (which is what glmIRLS uses.)
if isempty(pwts) && isscalar(nTrials)
    totalWeight = size(X,1);
elseif ~isempty(pwts) && isscalar(nTrials)
    totalWeight = sum(pwts);
elseif isempty(pwts) && ~isscalar(nTrials)
    totalWeight = sum(nTrials);
else
    totalWeight = sum(pwts .* nTrials);
end

% Convert lambda penalty to match total log-likelihood
lambda = lambda * totalWeight;

penalizedFitPartition = @(x,y,offset,pwts,n,wlsfit,b,active,mu,eta,sqrtvarFun) ...
    glmIRLSwrapper(x,y,distr,offset,pwts,dataClass,n, ...
    sqrtvarFun,linkFun,dlinkFun,ilinkFun,devFun,b,active,mu,muLims,wlsfit,nullDev,reltol);

penalizedFit = @(x,y,wlsfit,b,active,mu,eta) ...
    penalizedFitPartition(x,y,offset,pwts',nTrials,wlsfit,b,active,mu,eta,sqrtvarFun);

[B,Intercept,lambda,deviance] = ...
    lassoFit(X,Y,pwts,lambda,alpha,dfmax,standardize,reltol,lambdaMax*totalWeight,ever_active, ...
    penalizedFit,mu,eta,dataClass,userSuppliedLambda,nullDev,nullIntercept);

% Store the number of non-zero coefficients for each lambda.
df = sum(B~=0,1);

% The struct 'stats' will comprise the second return argument.
% Put place holders for ever-present fields to avoid undefined references
% and to secure the order we want in the struct.
stats = struct();
stats.Intercept      = [];
stats.Lambda         = [];
stats.Alpha          = alpha;
stats.DF             = [];
stats.Deviance       = [];
stats.PredictorNames = predictorNames;

% ---------------------------------------------------------
% If requested, use cross-validation to calculate
% Prediction Mean Squared Error for each lambda.
% ---------------------------------------------------------

if ~isequal(cvp,'resubstitution')
    % Replace dfmax with P, the number of predictors supplied at the
    % command line. dfmax might cause one fold to return empty values,
    % because no lambda satisfies the dfmax criteria, while other folds
    % return a numeric value. The lambda sequence has already been
    % pruned by dfmax, if appropriate, in the call to lassoFit above.
    cvfun = @(Xtrain,Ytrain,Xtest,Ytest) lassoFitAndPredict( ...
        Xtrain,Ytrain,Xtest,Ytest, ...
        lambda,alpha,P,standardize,reltol,ever_active, ...
        penalizedFitPartition,distr,link,linkFun,dlinkFun,sqrtvarFun, ...
        isCanonical, dlinkFunCanonical,devFun,dataClass);
    weights = pwts;
    if isempty(weights)
        weights = nan(size(X,1),1);
    end
    if isempty(offset) || isequal(offset,0)
        offset = nan(size(X,1),1);
    end
    if binomialTwoColumn
        response = [nTrials Y];
    else
        response = Y;
    end
    cvDeviance = crossval(cvfun,[weights(:) offset(:) X],response, ...
        'Partition',cvp,'Mcreps',mcreps,'Options',parallelOptions);
    % Scale single-fold deviances up to estimate whole data set deviance
    cvDeviance = bsxfun(@times,cvDeviance,repmat((size(X,1) ./ cvp.TestSize)', mcreps, 1));
    deviance = mean(cvDeviance);
    se = std(cvDeviance) / sqrt(size(cvDeviance,1));
    minDeviance = min(deviance);
    minIx = find(deviance == minDeviance,1);
    lambdaMin = lambda(minIx);
    minplus1 = deviance(minIx) + se(minIx);
    seIx = find((deviance(1:minIx) <= minplus1),1,'first');
    if isempty(seIx)
        lambdaSE = [];
    else
        lambdaSE = lambda(seIx);
    end
    
    % Deposit cross-validation results in struct for return value.
    stats.SE                = se;
    stats.LambdaMinDeviance = lambdaMin;
    stats.Lambda1SE         = lambdaSE;
    stats.IndexMinDeviance  = minIx;
    stats.Index1SE          = seIx;
    
    % Convert key lambda values determined by cross-validation
    % back to match average log-likelihood per observation.
    stats.LambdaMinDeviance = stats.LambdaMinDeviance / totalWeight;
    stats.Lambda1SE         = stats.Lambda1SE / totalWeight;
end

% ------------------------------------------
% Order results by ascending lambda
% ------------------------------------------

nLambda = length(lambda);
reverseIndices = nLambda:-1:1;
lambda = lambda(reverseIndices);
lambda = reshape(lambda,1,nLambda);
B = B(:,reverseIndices);
Intercept = Intercept(reverseIndices);
df = df(reverseIndices);
deviance = deviance(reverseIndices);
if ~isequal(cvp,'resubstitution')
    stats.SE               = stats.SE(reverseIndices);
    stats.IndexMinDeviance = nLambda - stats.IndexMinDeviance + 1;
    stats.Index1SE         = nLambda - stats.Index1SE + 1;
end

stats.Intercept = Intercept;
stats.Lambda = lambda;
stats.DF = df;
stats.Deviance = deviance;

% Convert lambda penalty back to match average deviance per observation.
stats.Lambda = stats.Lambda / totalWeight;

end %-main block

% ------------------------------------------
% SUBFUNCTIONS 
% ------------------------------------------

% ===================================================
%                  startingVals() 
% ===================================================

function mu = startingVals(distr,y,N)
% Find a starting value for the mean, avoiding boundary values
switch distr
case 'poisson'
    mu = y + 0.25;
case 'binomial'
    mu = (N .* y + 0.5) ./ (N + 1);
case {'gamma' 'inverse gaussian'}
    mu = max(y, eps(class(y))); % somewhat arbitrary
otherwise
    mu = y;
end
end %-startingVals

% ===================================================
%                diagnoseSeparation()
% ===================================================

function diagnoseSeparation(eta,y,N)
% Compute sample proportions, sorted by increasing fitted value
[x,idx] = sort(eta);
if ~isscalar(N)
    N = N(idx);
end
p = y(idx);
if all(p==p(1))   % all sample proportions are the same
    return
end
if x(1)==x(end)   % all fitted probabilities are the same
    return
end

noFront = 0<p(1) && p(1)<1;     % no "front" section as defined below
noEnd = 0<p(end) && p(end)<1;   % no "end" section as defined below
if p(1)==p(end) || (noFront && noEnd)
    % No potential for perfect separation if the ends match or neither
    % end is perfect
    return
end

% There is at least one observation potentially taking probability 0 or
% 1 at one end or the other with the data sorted by eta. We want to see
% if the data, sorted by eta (x) value, have this form:
%        x(1)<=...<=x(A)  <  x(A+1)=...=x(B-1)  <  x(B)<=...<=x(n)
% with   p(1)=...=p(A)=0                           p(B)=...=p(n)=1
% or     p(1)=...=p(A)=1                           p(B)=...=p(n)=0
%
% This includes the possibilities:
%     A+1=B  - no middle section
%     A=0    - no perfect fit at the front
%     B=n+1  - no perfect fit at the end
dx = 100*max(eps(x(1)),eps(x(end)));
n = length(p);
if noFront
    A = 0;
else
    A = find(p~=p(1),1,'first')-1;
    cutoff = x(A+1)-dx;
    A = sum(x(1:A)<cutoff);
end

if noEnd
    B = n+1;
else
    B = find(p~=p(end),1,'last')+1;
    cutoff = x(B-1)+dx;
    B = (n+1) - sum(x(B:end)>cutoff);
end

if A+1<B-1
    % There is a middle region with >1 point, see if x varies there
    if x(B-1)-x(A+1)>dx
        return
    end
end

% We have perfect separation that can be defined by some middle point
if A+1==B
    xmid = x(A) + 0.5*(x(B)-x(A));
else
    xmid = x(A+1);
    if isscalar(N)
        pmid = mean(p(A+1:B-1));
    else
        pmid = sum(p(A+1:B-1).*N(A+1:B-1)) / sum(N(A+1:B-1));
    end
end

% Create explanation part for the lower region, if any
if A>=1
    explanation = sprintf('\n   XB<%g: P=%g',xmid,p(1));
else
    explanation = '';
end

% Add explanation part for the middle region, if any
if A+1<B
    explanation = sprintf('%s\n   XB=%g: P=%g',explanation,xmid,pmid);
end
    
% Add explanation part for the upper region, if any
if B<=n
    explanation = sprintf('%s\n   XB>%g: P=%g',explanation,xmid,p(end));
end

warning(message('stats:lassoGlm:PerfectSeparation', explanation));
end %-diagnoseSeparation()

% ===============================================
%               glmProcessData() 
% ===============================================

function [x, y, offset, pwts, dataClass, N, binomialTwoColumn] = ...
    glmProcessData(x, y, distr, const, offset, pwts)

N = []; % needed only for binomial
binomialTwoColumn = false;

% Convert the two-column form of 'y', if supplied ('binomial' only).
if strcmp(distr,'binomial')
    if size(y,2) == 1
        % N will get set to 1 below
        if any(y < 0 | y > 1)
            error(message('stats:lassoGlm:BadDataBinomialFormat'));
        end
    elseif size(y,2) == 2
        binomialTwoColumn = true;
        y(y(:,2)==0,2) = NaN;
        N = y(:,2);
        y = y(:,1) ./ N;
        if any(y < 0 | y > 1)
            error(message('stats:lassoGlm:BadDataBinomialRange'));
        end
    else
        error(message('stats:lassoGlm:MatrixOrBernoulliRequired'));
    end
end

[anybad,~,y,x,offset,pwts,N] = dfswitchyard('statremovenan',y,x,offset,pwts,N);
if anybad > 0
    switch anybad
    case 2
        error(message('stats:lassoGlm:InputSizeMismatchX'))
    case 3
        error(message('stats:lassoGlm:InputSizeMismatchOffset'))
    case 4
        error(message('stats:lassoGlm:InputSizeMismatchPWTS'))
    end
end

% Extra screening for lassoglm (Infs and zero weights)
okrows = all(isfinite(x),2) & all(isfinite(y),2) & all(isfinite(offset));

if ~isempty(pwts)
    % This screen works on weights prior to stripping NaNs and Infs.
    if ~isvector(pwts) || ~isreal(pwts) || size(x,1) ~= length(pwts) || ...
            ~all(isfinite(pwts)) || any(pwts<0)
        error(message('stats:lassoGlm:InvalidObservationWeights'));
    end    
    okrows = okrows & pwts(:)>0;
    pwts = pwts(okrows);
end

% We need at least two observations after stripping NaNs and Infs and zero weights.
if sum(okrows)<2
    error(message('stats:lassoGlm:TooFewObservationsAfterNaNs'));
end

% Remove observations with Infs in the predictor or response
% or with zero observation weight.  NaNs were already gone.
x = x(okrows,:);
y = y(okrows);
if ~isempty(N) && ~isscalar(N)
    N = N(okrows);
end
if ~isempty(offset)
    offset = offset(okrows);
end

if isequal(const,'on')
    x = [ones(size(x,1),1) x];
end
dataClass = superiorfloat(x,y);
x = cast(x,dataClass);
y = cast(y,dataClass);

if isempty(offset), offset = 0; end
if isempty(N), N = 1; end

end %-glmProcessData()

% ===================================================
%             processLassoParameters() 
% ===================================================

function [X,Y,weights,nLambda,lambda,dfmax,cvp,mcreps,predictorNames,ever_active] = ...
    processLassoParameters(X,Y,weights, alpha, nLambda, lambdaRatio, lambda, dfmax, ...
    standardize, reltol, cvp, mcreps, predictorNames)

% === 'Weights' parameter ===
if ~isempty(weights)
        
    % Computations expect that weights is a row vector.
    weights = weights(:)';
    
end

[~,P] = size(X);

% If X has any constant columns, we want to exclude them from the
% coordinate descent calculations.  The corresponding coefficients
% will be returned as zero.
constantPredictors = (range(X)==0);
ever_active = ~constantPredictors;

% === 'Alpha' parameter ===

% Require 0 < alpha <= 1.
% "0" would correspond to ridge, "1" is lasso.
if ~isscalar(alpha) || ~isreal(alpha) || ~isfinite(alpha) || ...
        alpha <= 0 || alpha > 1
    error(message('stats:lassoGlm:InvalidAlpha'))
end

% === 'Standardize' option ===

% Require a logical value.
if ~isscalar(standardize) || (~islogical(standardize) && standardize~=0 && standardize~=1)
    error(message('stats:lassoGlm:InvalidStandardize'))
end

% === 'Lambda' sequence or 'NumLambda' and 'lambdaRatio' ===

if ~isempty(lambda)
    
    % Sanity check on user-supplied lambda sequence.  Should be non-neg real.
    if ~isreal(lambda) || any(lambda < 0)
        error(message('stats:lassoGlm:NegativeLambda'));
    end

    lambda = sort(lambda(:),1,'descend');
    
else
    
    % Sanity-check of 'NumLambda', should be positive integer.
    if ~isreal(nLambda) || ~isfinite(nLambda) || nLambda < 1
        error(message('stats:lassoGlm:InvalidNumLambda'));
    else
        nLambda = floor(nLambda);
    end
    
    % Sanity-checking of LambdaRatio, should be in [0,1).
    if ~isreal(lambdaRatio) || lambdaRatio <0 || lambdaRatio >= 1
        error(message('stats:lassoGlm:InvalidLambdaRatio'));
    end
end

% === 'RelTol' parameter ===
%
if ~isscalar(reltol) || ~isreal(reltol) || ~isfinite(reltol) || reltol <= 0 || reltol >= 1
    error(message('stats:lassoGlm:InvalidRelTol'));
end

% === 'DFmax' parameter ===
%
% DFmax is #non-zero coefficients 
% DFmax should map to an integer in [1,P] but we truncate if .gt. P
%
if isempty(dfmax)
    dfmax = P;
else
    if ~isscalar(dfmax)
        error(message('stats:lassoGlm:DFmaxBadType'));
    end
    try
        dfmax = uint32(dfmax);
    catch ME
        mm = message('stats:lassoGlm:DFmaxBadType');
        throwAsCaller(MException(mm.Identifier,'%s',getString(mm)));
    end
    if dfmax < 1
        error(message('stats:lassoGlm:DFmaxNotAnIndex'));
    else
        dfmax = min(dfmax,P);
    end
end

% === 'Mcreps' parameter ===
%
if ~isscalar(mcreps) || ~isreal(mcreps) || ~isfinite(mcreps) || mcreps < 1
    error(message('stats:lassoGlm:MCRepsBadType'));
end
mcreps = fix(mcreps);

% === 'CV' parameter ===
%

if isnumeric(cvp) && isscalar(cvp) && (cvp==round(cvp)) && (0<cvp)
    % cvp is a kfold value. It will be passed as such to crossval.
    if (cvp>size(X,1))
        error(message('stats:lassoGlm:InvalidCVforX'));
    end
    cvp = cvpartition(size(X,1),'Kfold',cvp);
elseif isa(cvp,'cvpartition')
    if strcmpi(cvp.Type,'resubstitution')
        cvp = 'resubstitution';
    elseif strcmpi(cvp.Type,'leaveout')
        error(message('stats:lassoGlm:InvalidCVtype'));
    elseif strcmpi(cvp.Type,'holdout') && mcreps<=1
        error(message('stats:lassoGlm:InvalidMCReps'));
    end
elseif strncmpi(cvp,'resubstitution',length(cvp))
    % This may have been set as the default, or may have been
    % provided at the command line.  In case it's the latter, we
    % expand abbreviations.
    cvp = 'resubstitution';
else
    error(message('stats:lassoGlm:InvalidCVtype'));
end
if strcmp(cvp,'resubstitution') && mcreps ~= 1
    error(message('stats:lassoGlm:InvalidMCReps'));
end

if isa(cvp,'cvpartition')
    if (cvp.N ~= size(X,1)) || (min(cvp.TrainSize) < 2)
        % We need partitions that match the total number of observations
        % (after stripping NaNs and Infs and zero observation weights), and
        % we need training sets with at least 2 usable observations.
        error(message('stats:lassoGlm:TooFewObservationsForCrossval'));
    end
end

% === 'PredictorNames' parameter ===
%
% If PredictorNames is not supplied or is supplied as empty, we just 
% leave it that way. Otherwise, confirm that it is a cell array of strings.
%
if ~isempty(predictorNames) 
    if ~iscellstr(predictorNames) || length(predictorNames(:)) ~= size(X,2)
        error(message('stats:lassoGlm:InvalidPredictorNames'));
    else
        predictorNames = predictorNames(:)';
    end
end

end %-processLassoParameters()

% ===================================================
%             glmProcessDistrAndLink()
% ===================================================

function [estdisp,sqrtvarFun,devFun,linkFun,dlinkFun,ilinkFun,link,mu,eta,muLims,...
    isCanonical,dlinkFunCanonical] = ...
    glmProcessDistrAndLink(y,distr,link,estdisp,N,dataClass)

switch distr
    case 'normal'
        canonicalLink = 'identity';
    case 'binomial'
        canonicalLink = 'logit';
    case 'poisson'
        canonicalLink = 'log';
    case 'gamma'
        canonicalLink = 'reciprocal';
    case 'inverse gaussian'
        canonicalLink = -2;
end

if isequal(link, 'canonical'), link = canonicalLink; end

switch distr
case 'normal'
    sqrtvarFun = @(mu) ones(size(mu));
    devFun = @(mu,y) (y - mu).^2;
    estdisp = 'on';
case 'binomial'
    sqrtN = sqrt(N);
    sqrtvarFun = @(mu) sqrt(mu).*sqrt(1-mu) ./ sqrtN;
    devFun = @(mu,y) 2*N.*(y.*log((y+(y==0))./mu) + (1-y).*log((1-y+(y==1))./(1-mu)));
case 'poisson'
    if any(y < 0)
        error(message('stats:lassoGlm:BadDataPoisson'));
    end
    sqrtvarFun = @(mu) sqrt(mu);
    devFun = @(mu,y) 2*(y .* (log((y+(y==0)) ./ mu)) - (y - mu));
case 'gamma'
    if any(y <= 0)
        error(message('stats:lassoGlm:BadDataGamma'));
    end
    sqrtvarFun = @(mu) mu;
    devFun = @(mu,y) 2*(-log(y ./ mu) + (y - mu) ./ mu);
    estdisp = 'on';
case 'inverse gaussian'
    if any(y <= 0)
        error(message('stats:lassoGlm:BadDataInvGamma'));
    end
    sqrtvarFun = @(mu) mu.^(3/2);
    devFun = @(mu,y) ((y - mu)./mu).^2 ./  y;
    estdisp = 'on';
otherwise
    error(message('stats:lassoGlm:BadDistribution'));
end

% Instantiate functions for one of the canned links, or validate a
% user-defined link specification.
[linkFun,dlinkFun,ilinkFun] = dfswitchyard('stattestlink',link,dataClass);

% Initialize mu and eta from y.
mu = startingVals(distr,y,N);
eta = linkFun(mu);

% Enforce limits on mu to guard against an inverse link that doesn't map into
% the support of the distribution.
switch distr
case 'binomial'
    % mu is a probability, so order one is the natural scale, and eps is a
    % reasonable lower limit on that scale (plus it's symmetric).
    muLims = [eps(dataClass) 1-eps(dataClass)];
case {'poisson' 'gamma' 'inverse gaussian'}
    % Here we don't know the natural scale for mu, so make the lower limit
    % small.  This choice keeps mu^4 from underflowing.  No upper limit.
    muLims = realmin(dataClass).^.25;
otherwise
    muLims = []; 
end

% These two quantities (isCanonical, dlinkFunCanonical) are not used by 
% glmfit but they are needed for calculation of lambdaMax in lassoglm.

isCanonical = isequal(link, canonicalLink);
[~, dlinkFunCanonical] = dfswitchyard('stattestlink', canonicalLink, dataClass);

end %-glmProcessDistrAndLink()

% ===================================================
%                      glmIRLS()
% ===================================================

function [b,mu,eta,varargout] = glmIRLS(x,y,distr,offset,pwts,dataClass,N, ...
    sqrtvarFun,linkFun,dlinkFun,ilinkFun,b,active,mu,muLims, ...
    wlsfit,nullDev,devFun,reltol)

wsIterationLimit = warning('off','stats:lassoGlm:IterationLimit');
wsPerfectSeparation = warning('off','stats:lassoGlm:PerfectSeparation');
wsBadScaling = warning('off','stats:lassoGlm:BadScaling');
cleanupIterationLimit = onCleanup(@() warning(wsIterationLimit));
cleanupPerfectSeparation = onCleanup(@() warning(wsPerfectSeparation));
cleanupBadScaling = onCleanup(@() warning(wsBadScaling));

if isempty(pwts)
    pwts = 1;
end

% Set up for iterations
iter = 0;
iterLim = 100;
warned = false;
seps = sqrt(eps);

% Match the convergence tolerance for the IRLS algorithm to the 
% command line parameter 'RelTol'.
convcrit = max(1e-6,2*reltol);

eta = linkFun(mu);

while iter <= iterLim
    iter = iter+1;

    % Compute adjusted dependent variable for least squares fit
    deta = dlinkFun(mu);
    z = eta + (y - mu) .* deta;

    % Compute IRLS weights the inverse of the variance function
    sqrtw = sqrt(pwts) ./ (abs(deta) .* sqrtvarFun(mu));
    
    % If the weights have an enormous range, we won't be able to do IRLS very
    % well.  The prior weights may be bad, or the fitted mu's may have too
    % wide a range, which is probably because the data do as well, or because
    % the link function is trying to go outside the distribution's support.
    wtol = max(sqrtw)*eps(dataClass)^(2/3);
    t = (sqrtw < wtol);
    if any(t)
        t = t & (sqrtw ~= 0);
        if any(t)
            sqrtw(t) = wtol;
            if ~warned
                warning(message('stats:lassoGlm:BadScaling'));
            end
            warned = true;
        end
    end

    b_old = b;
    [b,active] = wlsfit(z - offset, x, sqrtw.^2, b, active);

    % Form current linear predictor, including offset
    eta = offset + x * b;

    % Compute predicted mean using inverse link function
    mu = ilinkFun(eta);

    % Force mean in bounds, in case the link function is a wacky choice
    switch distr
    case 'binomial'
        if any(mu < muLims(1) | muLims(2) < mu)
        mu = max(min(mu,muLims(2)),muLims(1));
        end
    case {'poisson' 'gamma' 'inverse gaussian'}
        if any(mu < muLims(1))
        mu = max(mu,muLims(1));
        end
    end

    % Check stopping conditions
    % Convergence of the coefficients
    if (~any(abs(b-b_old) > convcrit * max(seps, abs(b_old)))) 
        break; 
    end
    % Proportion of explained deviance explained
    if sum(devFun(mu,y)) < (1.0e-3 * nullDev)
        break;
    end
     
end 

if iter > iterLim
    warning(message('stats:lassoGlm:IterationLimit'));
end

if iter>iterLim && isequal(distr,'binomial')
    diagnoseSeparation(eta,y,N);
end

varargout{1} = active;

end %-glmIRLS

% ===================================================
%                      glmIRLSwrapper()
% ===================================================

function [B,active,varargout] = glmIRLSwrapper(X,Y,distr,offset,pwts,dataClass,N, ...
    sqrtvarFun,linkFun,dlinkFun,ilinkFun,devFun,b,active,mu,muLims, ...
    wlsfit,nullDev,reltol)
% glmIRLSwrapper is a utility function for regularized GLM. It is called by 
% the regularization framework (lassoFit()), and adapts to the conventions
% and expectations of the IRLS implementation for GLM.
%
% Lasso always returns an intercept, but lassoFit passes in the predictor
% matrix 'X' without a column of ones. Unlike OLS, GLM cannot calve off 
% estimation of the predictor variable coefficients and calculate the
% intercept externally, post-fit.  This is because each IRLS step
% computes a local linear predictor (eta), and this requires an 
% intercept contribution. Therefore, prepend a column of ones before
% passing to glmIRLS. The variable B returned by glmIRLS will include
% the intercept term as the first element of B.  
X = [ones(size(X,1),1) X];

% glmIRLS assumes pwts=1 if no observation weights were supplied, 
% rather than empty, which is what lassoFit will feed in.
if isempty(pwts), pwts=1; end

[B,mu,eta,active] = glmIRLS(X,Y,distr,offset,pwts,dataClass,N, ...
    sqrtvarFun,linkFun,dlinkFun,ilinkFun,b,active,mu,muLims, ...
    wlsfit,nullDev,devFun,reltol);

deviance = sum(pwts.* devFun(mu,Y));


% Pull the intercept estimate out of B and return separately.
Intercept = B(1);
B = B(2:end);

extras.Intercept = Intercept;
extras.Deviance  = deviance;
varargout{1} = extras;
varargout{2} = mu;
varargout{3} = eta;

end %-glmIRLSwrapper

% ===================================================
%              lassoFitAndPredict() 
% ===================================================

function dev = lassoFitAndPredict(Xtrain,Ytrain,Xtest,Ytest, ...
    lambda,alpha,dfmax,standardize,reltol,ever_active, ...
    penalizedFitPartition,distr,link,linkFun,dlinkFun,sqrtvarFun, ...
    isCanonical, dlinkFunCanonical,devFun,dataClass)
%
% This function sets up the conditions for lasso fit and prediction from
% within crossvalidation.  It defines quantities as necessary for the
% current fold.

trainWeights = Xtrain(:,1);
% To conform with crossval syntax, weights are prepended to the predictor
% matrix. Empty weights are temporarily converted to a vector of NaNs.
% This clause restores them to empty.  Same occurs below for test weights.
if any(isnan(trainWeights))
    trainWeights = [];
end
trainOffset = Xtrain(:,2);
if any(isnan(trainOffset))
    trainOffset = 0;
end

Xtrain = Xtrain(:,3:end);
if size(Ytrain,2) == 2
    trainN = Ytrain(:,1);
    Ytrain = Ytrain(:,2);
else
    trainN = 1;
end

% These quantities depend on the particular set of responses used for the fit. 
% Within crossval, we work with a subset of the original data. Therefore, 
% these quantities must be redefined on each invocation.
mu = startingVals(distr,Ytrain,trainN);
eta = linkFun(mu);
if isequal(distr,'binomial')
    sqrtvarFun = @(mu) sqrt(mu).*sqrt(1-mu) ./ sqrt(trainN);
    devFun = @(mu,y) 2*trainN.*(y.*log((y+(y==0))./mu) + (1-y).*log((1-y+(y==1))./(1-mu)));
end

penalizedFit = @(x,y,wlsfit,b,active,mu,eta) penalizedFitPartition(x,y, ...
    trainOffset,trainWeights,trainN,wlsfit,b,active,mu,eta,sqrtvarFun);

[lambdaMax, nullDev, nullIntercept] = computeLambdaMax(Xtrain, Ytrain, trainWeights, ...
    alpha, standardize, distr, link, dlinkFun, trainOffset, isCanonical, dlinkFunCanonical, devFun);

% Will convert from lambda penalty matched to average log-likelihood
% per observation (which is what lasso uses) to an equivalent
% penalty for total log-likelihood (which is what glmIRLS uses.)
if isempty(trainWeights) && isscalar(trainN)
    totalWeight = size(Xtrain,1);
elseif ~isempty(trainWeights) && isscalar(trainN)
    totalWeight = sum(trainWeights);
elseif isempty(trainWeights) && ~isscalar(trainN)
    totalWeight = sum(trainN);
else
    totalWeight = sum(trainWeights .* trainN);
end

lambdaMax = lambdaMax * totalWeight;

[B,Intercept] = lassoFit(Xtrain,Ytrain, ...
    trainWeights,lambda,alpha,dfmax,standardize,reltol, ...
    lambdaMax,ever_active,penalizedFit,mu,eta,dataClass,true,nullDev,nullIntercept);
Bplus = [Intercept; B];

testWeights = Xtest(:,1);
if any(isnan(testWeights))
    testWeights = ones(size(Xtest,1),1);
end
testOffset = Xtest(:,2);
if any(isnan(testOffset))
    testOffset = 0;
end
Xtest = Xtest(:,3:end);
if size(Ytest,2) == 2
    testN = Ytest(:,1);
    Ytest = Ytest(:,2);
else
    testN = 1;
end

% For binomial with two-column form of response, the deviance function
% depends on the number of trials for each observation, which is provided 
% in testN. Within crossval, the test set is a subset of the original data.
% Therefore, the deviance function needs to be redefined for each test set.
if isequal(distr,'binomial')
    devFun = @(mu,y) 2*testN.*(y.*log((y+(y==0))./mu) + (1-y).*log((1-y+(y==1))./(1-mu)));
end

numFits = size(Bplus,2);
dev = zeros(1,numFits);
for i=1:numFits
    if ~isequal(testOffset,0)
        mu = glmval(Bplus(:,i), Xtest, link, 'Offset',testOffset);
    else
        mu = glmval(Bplus(:,i), Xtest, link);
    end
    di = devFun(mu,Ytest);
    dev(i) = sum(testWeights' * di);
end

end %-lassoFitAndPredict

% ===================================================
%                    lassoFit()                     
% ===================================================

function [B,Intercept,lambda,varargout] = ...
    lassoFit(X,Y,weights,lambda,alpha,dfmax,standardize,reltol, ...
    lambdaMax,ever_active,penalizedFit,mu,eta,dataClass,userSuppliedLambda,nullDev,nullIntercept)
%
% ------------------------------------------------------
% Perform model fit for each lambda and the given alpha
% ------------------------------------------------------

regressionType = 'GLM';

[~,P] = size(X);
nLambda = length(lambda);

% If X has any constant columns, we want to exclude them from the
% coordinate descent calculations.  The corresponding coefficients
% will be returned as zero.
constantPredictors = (range(X)==0);
ever_active = ever_active & ~constantPredictors;

% === weights and standardization ===
%
observationWeights = ~isempty(weights);
if ~isempty(weights)
    observationWeights = true;
    weights = weights(:)';
    % Normalize weights up front.
    weights = weights / sum(weights);
end

if standardize
    if ~observationWeights
        % Center and scale
        [X0,muX,sigmaX] = zscore(X,1);
        % Avoid divide by zero with constant predictors
        sigmaX(constantPredictors) = 1;
    else
        % Weighted center and scale
        muX = weights*X;
        X0 = bsxfun(@minus,X,muX);
        sigmaX = sqrt( weights*(X0.^2) );
        % Avoid divide by zero with constant predictors
        sigmaX(constantPredictors) = 1;
        X0 = bsxfun(@rdivide, X0, sigmaX);
    end
else
    switch regressionType
        case 'OLS'           
            if ~observationWeights
                % Center
                muX = mean(X,1);
                X0 = bsxfun(@minus,X,muX);
                sigmaX = 1;
            else
                % Weighted center
                muX = weights*X;
                X0 = bsxfun(@minus,X,muX);
                sigmaX = 1;
            end
        case 'GLM'
            X0 = X;
            % For GLM don't center until we get inside the IRLS
            sigmaX = 1;
            muX = zeros(1,size(X,2));
    end
end

% For OLS, Y is centered, for GLM it is not.
switch regressionType
    case 'OLS'
        if ~observationWeights
            muY = mean(Y);
        else
            muY = weights*Y;
        end
        Y0 = bsxfun(@minus,Y,muY);
    case 'GLM'
        Y0 = Y;
end

% Preallocate the returned matrix of coefficients, B, 
% and other variables sized to nLambda.

B = zeros(P,nLambda);

b = zeros(P,1,dataClass);

if nLambda > 0
    Extras(nLambda) = struct('Intercept', nullIntercept, 'Deviance', nullDev);
    for i=1:nLambda-1, Extras(i) = Extras(nLambda); end
    intercept = nullIntercept;
end

active = false(1,P);

for i = 1:nLambda
    
    lam = lambda(i);
    
    if lam >= lambdaMax
        continue;
    end
    
    % This anonymous function adapts the conventions of the IRLS algorithm
    % for glm fit to the coordinate descent algorithm for penalized WLS,
    % which is used within the IRLS algorithm. 
    wlsfit = @(x,y,weights,b,active) glmPenalizedWlsWrapper(y,x,b,active,weights,lam, ...
        alpha,reltol,ever_active);

    [b,active,extras,mu,eta] = penalizedFit(X0,Y0,wlsfit,[intercept;b],active,mu,eta);

    B(:,i) = b;
    
    Extras(i) = extras;
    
    % Halt if maximum model size ('DFmax') has been met or exceeded.
    if sum(active) > dfmax
        % truncate B and lambda output arguments
        lambda = lambda(1:(i-1));
        B = B(:,1:(i-1));
        Extras = Extras(:,1:(i-1));
        break
    end
    
    % Halt if we have exceeded a threshold on the percent of 
    % null deviance left unexplained.
    if ~(userSuppliedLambda || isempty(nullDev))
        if extras.Deviance < 1.0e-3 * nullDev
            lambda = lambda(1:i);
            B = B(:,1:i);
            Extras = Extras(:,1:i);
            break
        end
    end
    
end % of lambda sequence

% ------------------------------------------
% Unwind the centering and scaling (if any)
% ------------------------------------------

B = bsxfun(@rdivide, B, sigmaX');
B(~ever_active,:) = 0;

switch regressionType
    case 'OLS'
        Intercept = muY-muX*B;
    case 'GLM'
        Intercept = zeros(1,length(lambda));
        for i=1:length(lambda)
            Intercept(i) = Extras(i).Intercept;
        end
        if isempty(lambda)
            Intercept = [];
        else
            Intercept = Intercept - muX*B;
        end
end

% ------------------------------------------
% Calculate Mean Prediction Squared Error (or other GoF)
% ------------------------------------------

switch regressionType
    case 'OLS'
        Intercept = muY-muX*B;
        BwithI = [Intercept; B];
        fits = [ones(size(X,1),1) X]*BwithI;
        residuals = bsxfun(@minus, Y, fits);
        if ~observationWeights
            mspe = mean(residuals.^2);
        else
            % This line relies on the weights having been normalized.
            mspe = weights * (residuals.^2);
        end        
        varargout{1} = mspe;
    case 'GLM'
        deviance = zeros(1,length(lambda));
        for i=1:length(lambda)
            deviance(i) = Extras(i).Deviance;
        end
        if isempty(lambda)
            deviance = [];
        end
        varargout{1} = deviance;
end

end %-lassoFit

% ===================================================
%                 thresholdScreen() 
% ===================================================

function potentially_active = thresholdScreen(X0, wX0, Y0, ...
    b, active, threshold)
r = Y0 - X0(:,active)*b(active,:);
% We don't need the (b.*wX2)' term that one might expect, because it
% is zero for the inactive predictors.
potentially_active = abs(r' *wX0) > threshold;
end %-thresholdScreen

% ===================================================
%                 cdescentCycleNewCandidates() 
% ===================================================

function [b,active,wX2,wX2calculated,shrinkFactor] = ...
    cdescentCycleNewCandidates(X0, weights, wX0, wX2, wX2calculated, Y0, ...
    b, active, shrinkFactor, threshold, candidates)
% 
r = Y0 - X0(:,active)*b(active,:);
bold = b;

for j=find(candidates);
    % Regress j-th partial residuals on j-th predictor
    bj = wX0(:,j)' * r;
    
    margin = abs(bj) - threshold;
    
    % Soft thresholding
    if margin > 0
        if ~wX2calculated(j)
            wX2(j) = weights * X0(:,j).^2;
            wX2calculated(j) = true;
            shrinkFactor(j) = wX2(j) + shrinkFactor(j);
        end
        
        b(j) = sign(bj) .* margin ./ shrinkFactor(j);

        active(j) = true;
    end
    
    r = r - X0(:,j)*(b(j)-bold(j));
end
 
end %-cdescentCycleNewCandidates

% ===================================================
%                 cdescentCycleNoRecalc() 
% ===================================================

function [b,active] = ...
    cdescentCycleNoRecalc(X0, wX0, wX2, Y0, b, active, shrinkFactor, threshold)
% 
r = Y0 - X0(:,active)*b(active,:);
bwX2 = b.*wX2;
bold = b;

for j=find(active);
    % Regress j-th partial residuals on j-th predictor
    bj = wX0(:,j)' * r + bwX2(j);
    
    margin = abs(bj) - threshold;
    
    % Soft thresholding
    if margin > 0
        b(j) = sign(bj) .* margin ./ shrinkFactor(j);
    else
        b(j) = 0;
        active(j) = false;
    end
    
    r = r - X0(:,j)*(b(j)-bold(j));
end
 
end %-cdescentCycleNoRecalc

% ===================================================
%                    penalizedWls()                     
% ===================================================

function [b,varargout] = ...
    penalizedWls(X,Y,b,active,weights,lambda,alpha,reltol)

    weights = weights(:)';
    
    [~,P] = size(X);

    wX = bsxfun(@times,X,weights');
    
    wX2 = zeros(P,1);
    wX2(active) = (weights * X(:,active).^2)';
    wX2calculated = active;
    
    threshold = lambda * alpha;
    
    shrinkFactor = wX2 + lambda * (1-alpha);

    % Iterative coordinate descent until converged
    while true
        
        bold = b;
        old_active = active;

        [b,active] = cdescentCycleNoRecalc(X,wX,wX2,Y, b,active,shrinkFactor,threshold);

        if ~any( abs(b(old_active) - bold(old_active)) > reltol * max(1.0,abs(bold(old_active))) )
            % Cycling over the active set converged.
            % Do one full pass through the predictors.
            % If there is no predictor added to the active set, we're done.
            % Otherwise, resume the coordinate descent iterations.
            bold = b;
            potentially_active = thresholdScreen(X,wX,Y,b,active,threshold);
            new_candidates = potentially_active & ~active;
            if any(new_candidates)
                [b,new_active,wX2,wX2calculated,shrinkFactor] = ...
                    cdescentCycleNewCandidates(X,weights,wX,wX2,wX2calculated,Y, ...
                    b,active,shrinkFactor,threshold,new_candidates);
            else
                new_active = active;
            end

            if isequal(new_active, active)
                break
            else
                super_active = active | new_active;
                if ~any( abs(b(super_active) - bold(super_active)) > reltol * max(1.0,abs(bold(super_active))) )
                    % We didn't change the coefficients enough by the full pass
                    % through the predictors to warrant continuing the iterations.
                    % The active coefficients after this pass and the
                    % active coefficients prior to this pass differ.
                    % This implies that a coefficient changed between zero
                    % and a small numeric value.  There is no "fit-wise" reason
                    % to prefer the before and after coefficients, so we
                    % choose the more parsimonious alternative.
                    if sum(new_active) > sum(active)
                        b = bold;
                    else
                        active = new_active;
                    end
                    break
                else
                    active = new_active;
                end
            end
        end
        
    end
    
    varargout{1} = active;
        
end %-penalizedWls()
 
% ===================================================
%                    glmPenalizedWlsWrapper()                     
% ===================================================

function [b,varargout] = glmPenalizedWlsWrapper(X,Y,b,active,weights, ...
    lambda,alpha,reltol,ever_active)

% The coordinate descent in penalizedWls() assumes centered X and Y,
% and it assumes X does NOT have a column of ones.  

X0 = X(:,2:end);

weights = weights(:)';

normedWeights = weights / sum(weights);

% Center X
muX = normedWeights * X0;
X0 = bsxfun(@minus,X0,muX);

% Center Y.
muY = normedWeights * Y;
Y = Y - muY;

[bPredictors,varargout{1}] = penalizedWls(X0, Y, b(2:end), ...
    active,weights,lambda,alpha,reltol);

bPredictors(~ever_active,:) = 0;

% Since X and Y were centered for the WLS we can calculate the intercept
% after the fact.
Intercept = muY-muX*bPredictors;
b = [Intercept; bPredictors];

end %-glmPenalizedWlsWrapper()

function [lambdaMax, nullDev, nullIntercept] = computeLambdaMax(X, Y, weights, alpha, standardize, ...
    distr, link, dlinkFun, offset, isCanonical, dlinkFunCanonical, devFun)

% lambdaMax is the penalty term (lambda) beyond which coefficients
% are guaranteed to be all zero: there is no benefit to calculating
% penalized fits with lambda > lambdaMax.
%
% nullDev is the deviance of the fit using just a constant term.
%
% The input parameter 'devFun' is used only as a sanity check to see if glmfit 
% gets a plausible fit with intercept term only.

% Head off potential cruft in the command window.
wsIllConditioned2 = warning('off','stats:glmfit:IllConditioned');
wsIterationLimit = warning('off','stats:glmfit:IterationLimit');
wsPerfectSeparation = warning('off','stats:glmfit:PerfectSeparation');
wsBadScaling = warning('off','stats:glmfit:BadScaling');
cleanupIllConditioned2 = onCleanup(@() warning(wsIllConditioned2));
cleanupIterationLimit = onCleanup(@() warning(wsIterationLimit));
cleanupPerfectSeparation = onCleanup(@() warning(wsPerfectSeparation));
cleanupBadScaling = onCleanup(@() warning(wsBadScaling));

if ~isempty(weights)
    observationWeights = true;
    weights = weights(:)';        
    % Normalized weights are used for standardization and calculating lambdaMax.
    normalizedweights = weights / sum(weights);
else
    observationWeights = false;
end

[N,~] = size(X);

% If we were asked to standardize the predictors, do so here because
% the calculation of lambdaMax needs the predictors as we will use
% them to perform fits.

if standardize
    % If X has any constant columns, we want to protect against
    % divide-by-zero in normalizing variances.
    constantPredictors = (range(X)==0);

    if ~observationWeights
        % Center and scale
        [X0,~,~] = zscore(X,1);
    else
        % Weighted center and scale
        muX = normalizedweights * X;
        X0 = bsxfun(@minus,X,muX);
        sigmaX = sqrt( normalizedweights * (X0.^2) );
        % Avoid divide by zero with constant predictors
        sigmaX(constantPredictors) = 1;
        X0 = bsxfun(@rdivide, X0, sigmaX);
    end
else
    if ~observationWeights
        % Center
        muX = mean(X,1);
        X0 = bsxfun(@minus,X,muX);
    else
        % Weighted center
        muX = normalizedweights(:)' * X;
        X0 = bsxfun(@minus,X,muX);
    end
end

constantTerm = ones(length(Y),1);
if isscalar(offset)
    [coeffs,nullDev] = glmfit(constantTerm,Y,distr,'constant','off', ...
        'link',link, 'weights',weights);
    predictedMu = glmval(coeffs,constantTerm,link,'constant','off');
else
    [coeffs,nullDev] = glmfit(constantTerm,Y,distr,'constant','off', ...
        'link',link,'weights',weights,'offset',offset);
    predictedMu = glmval(coeffs,constantTerm,link,'constant','off','offset',offset);
end

nullIntercept = coeffs;

% Sanity check. With badly matched link / distr / data, glmfit may not 
% have been able to converge to a reasonble estimate.  If so, the inputs
% to lassoglm may constitute a difficult problem formulation with 
% unsatisfactory maximum likelihood solution.  Poor formulations
% have been observed with mismatched links (ie, 'reciprocal' link with
% the Poisson distribution, in place of canonical 'log').  As a screen for
% this contingency, calculate the deviance we would get using the scalar
% unmodeled mean for the response data. Call this quantity "muDev".  
% The value of muDev should be no better than the nullDev calculated
% by glmfit above (albeit it might be equal or nearly equal).
% If the value is better, warn that the computations are of uncertain validity.

if observationWeights
    muDev = weights * devFun(mean(Y)*ones(length(Y),1), Y);
else
    muDev = sum(devFun(mean(Y)*ones(length(Y),1), Y));
end
if (muDev - nullDev) / max([1.0 muDev nullDev]) < - 1.0e-4
    [~, lastid] = lastwarn;
    if strcmp(lastid,'stats:glmfit:BadScaling')
        % This reassignment of predicted values is not guaranteed to
        % improve matters, but is an attempt to generate a workable
        % sequence of lambda values. Note: Since it is a constant
        % value for all observations, observation weights are a wash.
        predictedMu = mean(Y)*ones(length(Y),1);
        warning(message('stats:lassoGlm:DifficultLikelihood'));
    end
end

if ~isCanonical
    X0 = bsxfun( @times, X0, dlinkFunCanonical(predictedMu) ./ dlinkFun(predictedMu) );
end

if ~observationWeights
    dotp = abs(X0' * (Y - predictedMu));
    lambdaMax = max(dotp) / (N*alpha);
else
    wX0 = bsxfun(@times, X0, normalizedweights');
    dotp = abs(sum(bsxfun(@times, wX0, (Y - predictedMu))));
    lambdaMax = max(dotp) / alpha;
end

end %-computeLambdaMax()

function lambda = computeLambdaSequence(lambdaMax, nLambda, lambdaRatio, LRdefault)

% Fill in the log-spaced sequence of lambda values.
        
if nLambda==1
    lambda = lambdaMax;
else
    % Fill in a number "nLambda" of smaller values, on a log scale.
    if lambdaRatio==0
        lambdaRatio = LRdefault;
        addZeroLambda = true;
    else
        addZeroLambda = false;
    end
    lambdaMin = lambdaMax * lambdaRatio;
    loghi = log(lambdaMax);
    loglo = log(lambdaMin);
    logrange = loghi - loglo;
    interval = -logrange/(nLambda-1);
    lambda = exp(loghi:interval:loglo)';
    if addZeroLambda
        lambda(end) = 0;
    else
        lambda(end) = lambdaMin;
    end
end

end %-computeLambdaSequence

