function [beta,r,J,Sigma,mse,errorModelInfo,robustw] = nlinfit(X,y,model,beta,varargin)
%NLINFIT Nonlinear least-squares regression.
%   BETA = NLINFIT(X,Y,MODELFUN,BETA0) estimates the coefficients of a
%   nonlinear regression function, using least squares estimation.  Y is a
%   vector of response (dependent variable) values.  Typically, X is a
%   design matrix of predictor (independent variable) values, with one row
%   for each value in Y and one column for each coefficient.  However, X
%   may be any array that MODELFUN is prepared to accept.  MODELFUN is a
%   function, specified using @, that accepts two arguments, a coefficient
%   vector and the array X, and returns a vector of fitted Y values.  BETA0
%   is a vector containing initial values for the coefficients.
%
%   [BETA,R,J,COVB,MSE] = NLINFIT(X,Y,MODELFUN,BETA0) returns the fitted
%   coefficients BETA, the residuals R, the Jacobian J of MODELFUN, the
%   estimated covariance matrix COVB for the fitted coefficients, and an
%   estimate MSE of the variance of the error term.  You can use these
%   outputs with NLPREDCI to produce confidence intervals for predictions,
%   and with NLPARCI to produce confidence intervals for the estimated
%   coefficients.  If you use a robust option (see below), you must use
%   COVB and may need MSE as input to NLPREDCI or NLPARCI to insure that
%   the confidence intervals take the robust fit properly into account.
%
%   [...] = NLINFIT(X,Y,MODELFUN,BETA0,OPTIONS) specifies control parameters
%   for the algorithm used in NLINFIT.  OPTIONS is a structure that can be
%   created by a call to STATSET.  Applicable STATSET parameters are:
%
%      'MaxIter'     - Maximum number of iterations allowed.  Defaults to 100.
%      'TolFun'      - Termination tolerance on the residual sum of squares.
%                      Defaults to 1e-8.
%      'TolX'        - Termination tolerance on the estimated coefficients
%                      BETA.  Defaults to 1e-8.
%      'Display'     - Level of display output during estimation.  Choices
%                      are 'off' (the default), 'iter', or 'final'.
%      'DerivStep'   - Relative difference used in finite difference gradient
%                      calculation.  May be a scalar, or the same size as
%                      the parameter vector BETA.  Defaults to EPS^(1/3).
%      'FunValCheck' - Check for invalid values, such as NaN or Inf, from
%                      the objective function.  'off' or 'on' (default).
%      'RobustWgtFun'- A weight function for robust fitting. Valid functions
%                      are 'bisquare', 'andrews', 'cauchy', 'fair', 'huber',
%                      'logistic', 'talwar', or 'welsch'. Default is '' (no
%                      robust fitting). Can also be a function handle that
%                      accepts a normalized residual as input and returns
%                      the robust weights as output.
%      'Robust'      - Robust will be removed in a future release. Use
%                      RobustWgtFun to invoke the robust fitting option.
%      'WgtFun'      - WgtFun will be removed in the future release. Use
%                      RobustWgtFun instead.
%      'Tune'        - The tuning constant used in robust fitting to normalize the
%                      residuals before applying the weight function.  A positive
%                      scalar.  The default value depends upon the weight function.
%                      This parameter is required if the weight function is
%                      specified as a function handle.
%
%   [BETA,R,J,COVB,MSE] = NLINFIT(X,Y,MODELFUN,BETA0,OPTIONS,...,'Weights',W)
%   nlinfit can accept an optional parameter name/value pair that specifies 
%   the observation weights:
% 
%     'Weights'        A vector of real positive weights the same size as Y, 
%                      each element of which specifies an observation weight. 
%                      Reducing the weight of an observation reduces the 
%                      influence of that observation on the fitted model. 
%                      This can also be specified as a function handle that 
%                      accepts a vector of predicted response values and 
%                      returns a vector of real positive weights as output. 
%                      Default is no weights.
% 
%    OPTIONS.RobustWgtFun must be [] when using observation weights. R and J
%    are weighted residuals and a weighted model function Jacobian respectively.
%    When J has full column rank, outputs MSE and R are related by: 
%                    MSE = (R'*R)/(n-p) 
%    where (n-p) = observations-parameters and COVB is related to J and MSE 
%    by:
%                   COVB = inv(J'*J) * MSE. 
%   
%   [BETA,R,J,COVB,MSE,ERRORMODELINFO] = NLINFIT(X,Y,MODELFUN,BETA0,OPTIONS,...
%                                            'ErrorModel',val1,'ErrorParameters',val2) 
%    nlinfit can also accept parameter name/value pairs that specify the error 
%    model and an initial estimate of error parameters to be used by nlinfit. 
%    These parameter name/value pairs are described below. 
% 
%     'ErrorModel'     A string specifying the form of the error term.
%                      Default is 'constant'. Each model defines the error
%                      using a standard zero mean and unit variance variable e
%                      with independent components, model function value f, 
%                      and one or two parameters a and b. The Choices for 
%                      val1 are:              
%                    
%                      'constant' (default)         y = f + a*e
%                      'proportional'               y = f + b*f*e
%                      'combined'                   y = f + (a+b*abs(f))*e
%
%                      If not specified, a 'constant' error model will be 
%                      used. The only allowed 'ErrorModel' when using 'Weights' 
%                      is 'constant'. OPTIONS.RobustWgtFun must be [] when 
%                      using error models other than 'constant'.
%
%   'ErrorParameters'  A numeric array containing initial estimates of the
%                      error model parameters of the chosen 'ErrorModel'. 
%                      Specify val2 as:   
%                     
%                      a       (default 1)       for 'constant'
%                      b       (default 1)       for 'proportional'
%                      [a b]   (default [1,1])   for 'combined'
% 
%                      For example, if 'ErrorModel' is 'combined', one could 
%                      use 'ErrorParameters' as [1,2]. If not specified, the 
%                      default values mentioned in brackets above will be used 
%                      as initial estimates.
%
%   The output ERRORMODELINFO is a structure with the following fields:
%
%                      ErrorModel   Chosen error model (default = 'constant')
%                 ErrorParameters   Estimated error parameters
%                   ErrorVariance   A function handle that accepts a n by p 
%                                   matrix X and computes a n by 1 vector of 
%                                   error variances using the specified error 
%                                   model.
%                             MSE   Mean squared error.
%                  ScheffeSimPred   Scheffe parameter for a simultaneous 
%                                   prediction interval when using this error 
%                                   model. 
%                  WeightFunction   True if a custom weight function was used
%                                   previously in nlinfit.
%                    FixedWeights   True if fixed weights were used previously
%                                   in nlinfit.
%            RobustWeightFunction   True if a robust fitting was used previously
%                                   in nlinfit.
%
%   When 'ErrorModel' is 'combined' or 'proportional', R and J can no
%   longer be interpreted as model fit residuals and model function
%   Jacobian respectively. When J has full column rank, outputs MSE and R
%   are related by:
%                    MSE = (R'*R)/(n-p) where (n-p) = observations-parameters 
%   and COVB is related to J and MSE by:
%                   COVB = inv(J'*J) * MSE 
%   
%   NLINFIT treats NaNs in Y or MODELFUN(BETA0,X) as missing data, and
%   ignores the corresponding observations.
%
%   Examples:
%
%      MODELFUN can be an anonymous function:
%         load reaction;
%         fun = @(beta,x)(beta(1)*x(:,2) - x(:,3)/beta(5)) ./ ...
%                        (1+beta(2)*x(:,1)+beta(3)*x(:,2)+beta(4)*x(:,3));
%         beta = nlinfit(reactants,rate,fun,beta);
%
%      MODELFUN can be a function on the path specified using @:
%         load reaction;
%         beta = nlinfit(reactants,rate,@hougen,beta);
%
%      For an example of weighted fitting, see the example "Weighted
%      Nonlinear Regression".
%
%   See also NLPARCI, NLPREDCI, NLMEFIT, NLINTOOL, STATSET.

%   References:
%      [1] Seber, G.A.F, and Wild, C.J. (1989) Nonlinear Regression, Wiley.

%   NLINFIT can be used to make a weighted fit with known weights:
%
%      load reaction;
%      w = [8 2 1 6 12 9 12 10 10 12 2 10 8]'; % some example known weights
%      ratew = sqrt(w).*rate;
%      mymodelw = @(beta,X) sqrt(w).*mymodel(beta,X);
%
%      [betaw,residw,Jw] = nlinfit(reactants,ratew,mymodelw,beta);
%      betaciw = nlparci(betaw,residw,Jw);
%      [ratefitw, deltafitw] = nlpredci(@mymodel,reactants,betaw,residw,Jw);
%      rmse = norm(residw) / (length(w)-length(rate))
%
%   Predict at the observed x values.  However, the prediction band
%   assumes a weight (measurement precision) of 1 at these points.
%
%      [ratepredw, deltapredw] = ...
%            nlpredci(@mymodel,reactants,betaw,residw,Jw,[],[],'observation');

%   Copyright 1993-2014 The MathWorks, Inc.


if nargin < 4
    error(message('stats:nlinfit:TooFewInputs'));
elseif ~isvector(y)
    error(message('stats:nlinfit:NonVectorY'));
end

% Parse input arguments
[errormodel, weights, errorparam, options, iterative, maxweight] = parseInVarargin(varargin(:));

% Check sizes of the model function's outputs while initializing the fitted
% values, residuals, and SSE at the given starting coefficient values.
model = fcnchk(model);
try
    yfit = model(beta,X);
catch ME
    if isa(model, 'inline')
        m = message('stats:nlinfit:ModelInlineError');
        throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
    elseif strcmp('MATLAB:UndefinedFunction', ME.identifier) ...
            && ~isempty(strfind(ME.message, func2str(model)))
        error(message('stats:nlinfit:ModelFunctionNotFound', func2str( model )));
    else
        m = message('stats:nlinfit:ModelFunctionError',func2str(model));
        throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
    end
end
if ~isequal(size(yfit), size(y))
    expSizeStr = sprintf('%-d-by-',size(y));
    actSizeStr = sprintf('%-d-by-',size(yfit));
    error(message('stats:nlinfit:WrongSizeFunOutput',expSizeStr(1:end-4),actSizeStr(1:end-4)));
end

% Check weights
nanweights = checkWeights(weights, yfit, y);

if strcmp(errormodel, 'exponential')
    % Transform model.
    [y, model] = applyLogTransformation(y, model);
    % Transform yfit as well.
    [yfit, ~] = applyLogTransformation(yfit, []);
end

% Find NaNs in either the responses or in the fitted values at the starting
% point.  Since X is allowed to be anything, we can't just check for rows
% with NaNs, so checking yhat is the appropriate thing.  Those positions in
% the fit will be ignored as missing values.  NaNs that show up anywhere
% else during iteration will be treated as bad values.
nans = (isnan(y(:)) | isnan(yfit(:)) | nanweights(:)); % a col vector
r = y(:) - yfit(:);
r(nans) = [];
n = numel(r);
p = numel(beta);
sse = r'*r;

% After removing NaNs in either the responses or in the fitted values at 
% the starting point, if n = 0 then stop execution.
if ( n == 0 )
    error(message('stats:nlinfit:NoUsableObservations'));
end

funValCheck = strcmp(options.FunValCheck, 'on');
if funValCheck && ~isfinite(sse), checkFunVals(r); end

% Set the level of display
switch options.Display
    case 'off',    verbose = 0;
    case 'notify', verbose = 1;
    case 'final',  verbose = 2;
    case 'iter',   verbose = 3;
end
maxiter = options.MaxIter;
robustw = [];

if isempty(options.RobustWgtFun)
    % display format for non-robust fit
    if verbose > 2 % iter
        disp(' ');
        disp('                                     Norm of         Norm of');
        disp('   Iteration             SSE        Gradient           Step ');
        disp('  -----------------------------------------------------------');
        disp(sprintf('      %6d    %12g',0,sse)); %#ok<DSPS>
    end
    
    if iterative
        % Iteratively re-weighted fitting
        [beta,J,cause,errorparam,fullr] = nlweightedfit(X,y,beta,model,options,verbose,maxiter, errormodel, errorparam, weights, maxweight);
    else
        % Known weights                
        % Cap the weights at max weight
        weights(weights > maxweight) = maxweight;
        
        w = sqrt(weights);
        yw = y .* w;
        modelw = @(b,x) w.*model(b,x);
        [beta,J,~,cause,fullr] = LMfit(X,yw, modelw,beta,options,verbose,maxiter);
    end
    
else
    % Allow for empty or scalar weights
    if isempty(weights)
        weights = ones(size(y));
    elseif isscalar(weights)
        weights = repmat(weights,size(y));
    end

    % Do a preliminary fit just to get residuals and leverage from the
    % least squares coefficient.  We won't count this against the iteration
    % limit.
    [beta_ls,J] = LMfit(X,y, model,beta,options,0,maxiter,weights);
    res = y - model(beta_ls,X);
    res(isnan(res)) = [];
    ols_s = norm(res) / sqrt(max(1,length(res)-numel(beta)));
    
    % display format for robust fit
    % Please note there are two loops for robust fit. It would be very
    % confusing if we display all iteration results. Instead, only the last
    % step of each inner loop (LM fit) will be output.
    if verbose > 2 % iter
        disp(' ');
        disp(getString(message('stats:nlinfit:IterationsRecalculateTheRobustWeights')));
        disp(' ');
        disp('   Iteration             SSE ');
        disp('  -----------------------------');
        disp(sprintf('      %6d    %12g',0,sse)); %#ok<DSPS>
    end
    [beta,J,sig,cause,fullr,robustw] = nlrobustfit(X,y,beta,model,J,ols_s,options,verbose,maxiter,weights);
end;

switch(cause)
    case 'maxiter'
        warning(message('stats:nlinfit:IterationLimitExceeded'));
    case 'tolx'
        if verbose > 1 % 'final' or 'iter'
            disp(getString(message('stats:nlinfit:TerminatedRelativeNorm')));
        end
    case 'tolfun'
        if verbose > 1 % 'final' or 'iter'
            disp(getString(message('stats:nlinfit:TerminatedRelativeChangeInSSE')));
        end
    case 'stall'
        warning(message('stats:nlinfit:UnableToDecreaseSSE'));
end

% If the Jacobian is ill-conditioned, then it's likely that two parameters
% are aliased and the estimates will be highly correlated.  Prediction at
% new x values not in the same column space is dubious.  NLPARCI will have
% trouble computing CIs because the inverse of J'*J is difficult to get
% accurately. NLPREDCI will have the same difficulty, and in addition,
% will in effect end up taking the difference of two very large, but nearly
% equal, variance and covariance terms, lose precision, and so the
% prediction bands will be erratic. It may also be that the Jacobian has
% one or more columns of zeros, meaning model is constant with respect to
% one or more parameters.  This may be because those parameters are not
% even in the expression in the model function, or they are multiplied by
% another param that is estimated at exactly zero (or something similar),
% or because some part of the model function is underflowing, making it a
% constant zero.

% Final J might have Inf/NaN values or may not be of the right size due to
% NaN removal in intermediate steps if FunValCheck was 'off'. If so, this
% is an error since no meaningful statistics can be computed from such a J.
% If FunValCheck is 'on', we know that the Jacobian must be good.
if ( funValCheck == false )
    nonfiniteJ = any(~isfinite(J(:)));
    wrongsizeJ = ( size(J,1) ~= n || size(J,2) ~= p );    
    if ( nonfiniteJ || wrongsizeJ )
        warning(message('stats:nlinfit:BadFinalJacobian'));         
        % [beta,r,J,Sigma,mse,errorModelInfo,robustw] = nlinfit(...)
        % beta will *not* have the converged value but some intermediate
        % value. Get "full" r. Fill in NaN's for J, Sigma and mse. 
        % Use sch = p + 1 to get errorModelInfo. Make robustw all NaN's.
        r = fullr;
        J = NaN(numel(fullr),p);        
        Sigma = NaN(p,p);
        mse = NaN;
        sch = p + 1; % conservative setting.
        errorModelInfo = fillErrorModelInfo(sch, errormodel, errorparam, mse, model, beta, weights, options);
        robustw = NaN(size(robustw));
        return;
    end
end

% Get QR factorization of J and mark if J is ill-conditioned.
[Q,R] = qr(J,0);
illConditionedJ = false;
if n <= p
    illConditionedJ = true;
    warning(message('stats:nlinfit:Overparameterized'));
elseif condest(R) > 1/sqrt(eps(class(beta)))
    illConditionedJ = true;
    if any(all(abs(J)<sqrt(eps(norm(J,1))),1),2) % one or more columns of zeros
        warning(message('stats:nlinfit:ModelConstantWRTParam'));
    else
        % no columns of zeros
        warning(message('stats:nlinfit:IllConditionedJacobian'));
    end
end

% We have beta, J, fullr and we need [beta,r,J,Sigma,mse,errorModelInfo].
% Inspect nargout and compute required quantities.
if nargout > 1
    % Return residuals and Jacobian that have missing values where needed.
    r = fullr;
    JJ(~nans,:) = J;
    JJ(nans,:) = NaN;
    J = JJ;
end

if nargout > 3   
    % Get rankJ, pinvJTJ and VQ using either:
    % (a) SVD - if J is ill-conditioned or
    % (b)  QR - if J is well-conditioned. In this case, we can set VQ equal 
    %           to Q from QR factorization of J to test the Scheffe parameter 
    %           condition.
    if illConditionedJ
        % Call isEstimable.
        TolSVD = eps(class(beta));
        [~,rankJ,pinvJTJ,~,~,VQ] = internal.stats.isEstimable(eye(numel(beta)),'DesignMatrix',J(~nans,:),'TolSVD',TolSVD);
    else
        % Use existing QR factorization.
        rankJ = size(J,2);
        Rinv = R \ eye(size(R));
        pinvJTJ = Rinv*Rinv';
        VQ = Q;
    end
    
    % Get MSE using either residuals or from Robust fit.
    if isempty(options.RobustWgtFun)
        % Can use residuals to get MSE. 
        % The denominator degrees of freedom is (n-rankJ).
        mse = sum(abs(r(~nans)).^2)/(n-rankJ);
    else
        % Using Robust fitting.
        mse = sig.^2;
    end
    
    % Get Sigma, the co-variance matrix of beta. 
    % Sigma may be singular when J is ill-conditioned.
    Sigma = mse*pinvJTJ;    
end

if nargout > 5    
    % Estimate the Scheffe parameter for simultaneous prediction intervals.
    % We will pass in VQ from above to avoid an SVD computation in getscheffeparam.
    if ~isa(weights,'function_handle') && isnumeric(weights) && ~isscalar(weights) && ~all(weights==1)
        % Fixed weight vector. Conservative setting is: sch = rankJ + 1. See if we 
        % can do better than this.
        EVFit = 1./weights(~nans); EVPred = ones(size(EVFit));
        sch = internal.stats.getscheffeparam('WeightedJacobian',J(~nans,:),'ErrorVarianceFit',EVFit,'ErrorVariancePred',EVPred,'Intopt','observation','VQ',VQ);         
    else
       % Either weights is a function_handle or errormodel = 'constant','proportional','combined' or 'exponential'. 
        sch = internal.stats.getscheffeparam('WeightedJacobian',J(~nans,:),'Intopt','observation','VQ',VQ); 
    end   
    
    % Create the structure ErrorModelInfo.
    errorModelInfo = fillErrorModelInfo(sch, errormodel, errorparam, mse, model, beta, weights, options);        
end

end % function nlinfit

%==== subfunction fillErrorModelInfo ====
function S = fillErrorModelInfo(sch, errormodel, errorparam, mse, model, beta, weights, options)
    %  Initialize output structure with the following fields:
    %  ErrorModel           - Name of the error model.
    %  ErrorParameters      - Parameters of the error model.
    %  ErrorVariance        - Variance function for the error model.
    %  MSE                  - Mean squared error.  
    %  ScheffeSimPred       - Scheffe parameter for simultaneous prediction intervals.
    %  WeightFunction       - True if using 'Weights' with a weight function.
    %  FixedWeights         - True if using 'Weights' with a fixed weight vector.
    %  RobustWeightFunction - True if using Robust fitting with options.RobustWgtFun.
    S = struct('ErrorModel',[],'ErrorParameters',[],'ErrorVariance',[],'MSE',[],'ScheffeSimPred',sch,'WeightFunction',false,'FixedWeights',false,'RobustWeightFunction',false);
     
    % Set the fields of S. We will burn in the values of beta, model, mse, 
    % errorparam and weights when creating function handle ErrorVariance.
    switch lower(errormodel)        
        case 'combined'
            % Note that errorparam is already set during model fitting.
            S.ErrorModel = 'combined';
            S.ErrorParameters = errorparam;
            S.MSE = mse;
            S.ErrorVariance = @(x) mse * ( errorparam(1) + errorparam(2)*abs(model(beta,x)) ).^2;
            
        case 'proportional'
            % We must set errorparam explicitly.
            errorparam = sqrt(mse);
            S.ErrorModel = 'proportional';
            S.ErrorParameters = errorparam;
            S.MSE = mse;
            S.ErrorVariance = @(x) mse * ( abs(model(beta,x)) ).^2;
            
        case 'exponential'
            % This is the 'constant' variance model but *after* log transformation.
            errorparam = sqrt(mse);
            S.ErrorModel = 'exponential';
            S.ErrorParameters = errorparam;
            S.MSE = mse;
            S.ErrorVariance = @(x) mse * ones(size(x,1),1);
            
        case 'constant'            
            if isa(weights, 'function_handle')
                % If weights is a function handle, we need to compute 
                % sigma^2 * diag(W^{-1}) as the ErrorVariance function.
                errorparam = sqrt(mse);
                S.ErrorModel = 'constant';
                S.ErrorParameters = errorparam;
                S.MSE = mse;
                S.ErrorVariance = @(x) mse * (1./weights(model(beta,x)));
                S.WeightFunction = true;
            elseif isnumeric(weights) && ~isscalar(weights) && ~all(weights==1)
                % Known weight vector.
                errorparam = sqrt(mse);
                S.ErrorModel = 'constant';
                S.ErrorParameters = errorparam;
                S.MSE = mse;
                S.ErrorVariance = @(x) mse * ones(size(x,1),1);
                S.FixedWeights = true;
            else
                % The 'constant' variance model, maybe using robust regression.
                errorparam = sqrt(mse);
                S.ErrorModel = 'constant';
                S.ErrorParameters = errorparam;
                S.MSE = mse;
                S.ErrorVariance = @(x) mse * ones(size(x,1),1);                
            end
            if ~isempty(options.RobustWgtFun)
                % Robust regression.
                S.RobustWeightFunction = true;
            end
        otherwise
                % This error model is not known - we should never get here.
                error(message('stats:nlinfit:InvalidErrorModel'));
    end
    
end % End of fillErrorModelInfo.

%----------------------------------------------------------------------
function  [beta,J,iter,cause,fullr] = LMfit(X,y, model,beta,options,verbose,maxiter,weights)
% Levenberg-Marquardt algorithm for nonlinear regression

% Set up convergence tolerances from options.
betatol = options.TolX;
rtol = options.TolFun;
fdiffstep = options.DerivStep;
if isscalar(fdiffstep)
    fdiffstep = repmat(fdiffstep, size(beta));
else
    % statset ensures fdiffstep is not a matrix.
    % Here, we ensure fdiffstep has the same shape as beta.
    fdiffstep = reshape(fdiffstep, size(beta));
end
funValCheck = strcmp(options.FunValCheck, 'on');

% Set initial weight for LM algorithm.
lambda = .01;

% Set the iteration step
sqrteps = sqrt(eps(class(beta)));

p = numel(beta);

if nargin<8 || isempty(weights)
    sweights = ones(size(y));
else
    sweights = sqrt(weights);
end

% treatment for nans
yfit = model(beta,X);
fullr = sweights(:) .* (y(:) - yfit(:));
nans = isnan(fullr); % a col vector
r = fullr(~nans);
sse = r'*r;

zerosp = zeros(p,1,class(r));
iter = 0;
breakOut = false;
cause = '';

while iter < maxiter
    iter = iter + 1;
    betaold = beta;
    sseold = sse;
    
    % Compute a finite difference approximation to the Jacobian
    J = getjacobian(beta,fdiffstep,model,X,yfit,nans,sweights);
    
    %% AR filtering step
    if(iter>maxiter/2)
        [ar, p] = nirs.math.robust_ar_fit(r,10);
        f = [1; -ar(2:end)];
    else
        f=1;
    end
    
    J=filter(f,1,J);
    r=filter(f,1,r);
    
    % Levenberg-Marquardt step: inv(J'*J+lambda*D)*J'*r
    diagJtJ = sum(abs(J).^2, 1);
    if funValCheck && ~all(isfinite(diagJtJ)), checkFunVals(J(:)); end
    Jplus = [J; diag(sqrt(lambda*diagJtJ))];  % AR-filtered version
    rplus = [r; zerosp];
    step = pinv(Jplus)* rplus;
   % step=robustfit(Jplus,rplus,'bisquare',4.68,'off',ones(size(rplus)),0);
    
    beta(:) = beta(:) + step;
    
    % Evaluate the fitted values at the new coefficients and
    % compute the residuals and the SSE.
    yfit = model(beta,X);
    fullr = sweights(:) .* (y(:) - yfit(:));
    r = fullr(~nans);
    sse = r'*r;
    if funValCheck && ~isfinite(sse), checkFunVals(r); end
    % If the LM step decreased the SSE, decrease lambda to downweight the
    % steepest descent direction.  Prevent underflowing to zero after many
    % successful steps; smaller than eps is effectively zero anyway.
    if sse < sseold
        lambda = max(0.1*lambda,eps);
        
        % If the LM step increased the SSE, repeatedly increase lambda to
        % upweight the steepest descent direction and decrease the step size
        % until we get a step that does decrease SSE.
    else
        while sse > sseold
            lambda = 10*lambda;
            if lambda > 1e16
                breakOut = true;
                break
            end
            Jplus = [J; diag(sqrt(lambda*sum(J.^2,1)))];
            step = pinv(Jplus) * rplus;
            beta(:) = betaold(:) + step;
            yfit = model(beta,X);
            fullr = sweights(:) .* (y(:) - yfit(:));
            r = fullr(~nans);
            sse = r'*r;
            if funValCheck && ~isfinite(sse), checkFunVals(r); end
        end
    end
    if verbose > 2 % iter
        disp(sprintf('      %6d    %12g    %12g    %12g', ...
            iter,sse,norm(2*r'*J),norm(step))); %#ok<DSPS>
    end
    
    % Check step size and change in SSE for convergence.
    if norm(step) < betatol*(sqrteps+norm(beta))
        cause = 'tolx';
        break
    elseif abs(sse-sseold) <= rtol*sse
        cause = 'tolfun';
        break
    elseif breakOut
        cause = 'stall';
        break
    end
end
if (iter >= maxiter)
    cause = 'maxiter';
end
end % function LMfit

%--------------------------------------------------------------------------
function checkFunVals(v)
% check if the function has finite output
if any(~isfinite(v))
    error(message('stats:nlinfit:NonFiniteFunOutput'));
end
end % function checkFunVals

%--------------------------------------------------------------------------
function [beta,Jw,cause,errorModelParam,fullr]=nlweightedfit(x,y,beta,model,options,verbose,maxiter,errorModel,errorModelParam,wgtfun,maxweight)
% iteratively re-weighted fit

betatol = options.TolX;
sqrteps = sqrt(eps(class(beta)));

if strcmpi(errorModel, 'combined')
    % Do a prelimary unweighted fit before calculating [a b] for combined error model
    beta = LMfit(x,y,model,beta,options,0,maxiter);
end

yfit = model(beta,x);
fullr = y(:) - yfit(:);
ok = ~isnan(fullr);

% Main loop of repeated nonlinear fits, adjust weights each time
totiter = 0;
w = NaN(size(y));
while maxiter>0
    beta0=beta;
    
    % Update weights
    yfit = model(beta,x);
    switch errorModel
        case 'proportional'
            w = 1./abs(yfit);
        case 'combined'
            ab = abs(fminsearch(@(ab) error_ab(ab,y,yfit),errorModelParam));
            w = 1./abs(ab(1)+ab(2)*abs(yfit));
            errorModelParam = ab;
        case 'constant'
            % You dont get here until weights is a function handle
            w = realsqrt(wgtfun(yfit));
    end
    
    % Cap the weights at max weight
    w(w > sqrt(maxweight)) = sqrt(maxweight);
    
    % this is weighted nlinfit
    yw = y .* w;
    modelw = @(b,x) w.*model(b,x);
    [beta,Jw,lsiter,cause,fullr] = LMfit(x,yw,modelw,beta0,options,0,maxiter); % 6th arg always silences display
    totiter = totiter + lsiter;
    maxiter = maxiter - lsiter;
    r = fullr(ok);
    
    % if there is no change in any coefficient, the iterations stop.
    if norm(beta-beta0) < betatol*(sqrteps+norm(beta))
        cause = 'tolx';
        break;
    end
    
    if verbose > 2 % iter
        disp(sprintf('      %6d    %12g', ...
            totiter, r'*r)); %#ok<DSPS>
    end
end

% this is a warning about the non-convergence
if maxiter<=0
    cause = 'maxiter';
end

end % function nlweighted fit

%--------------------------------------------------------------------------
function [beta,J,sig,cause,fullr,w]=nlrobustfit(x,y,beta,model,J,ols_s,options,verbose,maxiter,pweights)
% nonlinear robust fit

tune = options.Tune;
WgtFun = options.RobustWgtFun;
[WgtFun,tune] = nirs.math.NonLinearModel.statrobustwfun(WgtFun,tune);

yfit = model(beta,x);
fullr = y - yfit;
ok = ~isnan(fullr(:));
r = fullr(:);
pweights = pweights(:);
r = sqrt(pweights(ok)).*r(ok);
Delta = sqrt(eps(class(beta)));

% Adjust residuals using leverage, as advised by DuMouchel & O'Brien
% Compute leverage based on X, the Jacobian
[Q,~]=qr(J,0);
h = min(.9999, sum(Q.*Q,2));

% Compute adjustment factor
adjfactor = 1 ./ sqrt(1-h);

radj = r .* adjfactor;

% If we get a perfect or near perfect fit, the whole idea of finding
% outliers by comparing them to the residual standard deviation becomes
% difficult.  We'll deal with that by never allowing our estimate of the
% standard deviation of the error term to get below a value that is a small
% fraction of the standard deviation of the raw response values.
tiny_s = 1e-6 * std(y);
if tiny_s==0
    tiny_s = 1;
end

% Main loop of repeated nonlinear fits, adjust weights each time
totiter = 0;
w = NaN(size(y));
while maxiter>0
    beta0=beta;
    s = madsigma(radj, numel(beta)); % robust estimate of sigma for residual
    
    % Compute robust weights based on current residuals
    w(ok) = feval(WgtFun, radj/(max(s,tiny_s)*tune));
    
    % this is the weighted nlinfit
    allweights = w(:).*pweights;
    allweights = reshape(allweights,size(y));
    [beta,~,lsiter,cause] = LMfit(x,y,model,beta0,options,0,maxiter,allweights); % 6th arg always silences display
    totiter = totiter + lsiter;
    maxiter = maxiter - lsiter;
    yfit = model(beta,x);
    fullr = y - yfit;
    r = fullr(:);
    r = sqrt(pweights(ok)).*r(ok);
    radj = r .* adjfactor;
    
    % if there is no change in any coefficient, the iterations stop.
    if  all(abs(beta-beta0) < Delta*max(abs(beta),abs(beta0)))
        break;
    end
    
    if verbose > 2 % iter
        disp(sprintf('      %6d    %12g', ...
            totiter, r'*r)); %#ok<DSPS>
    end
end

% this is a warning about the non-convergence
if maxiter<=0
    cause = 'maxiter';
end

% We need the Jacobian at the final coefficient estimates, but not the J1
% version returned by LMfit because it has robust weights included
fdiffstep = options.DerivStep;
if isscalar(fdiffstep)
    fdiffstep = repmat(fdiffstep,size(beta));
end
J = getjacobian(beta,fdiffstep,model,x,yfit,~ok,reshape(sqrt(pweights),size(yfit)));

% Compute MAD of adjusted residuals after dropping p-1 closest to 0
p = numel(beta);
n = length(radj);
mad_s = madsigma(radj, p);

% Compute a robust scale estimate for the covariance matrix
sig = nirs.math.NonLinearModel.statrobustsigma(WgtFun,radj,p,mad_s,tune,h);

% Be conservative by not allowing this to be much bigger than the ols value
% if the sample size is not large compared to p^2
sig = max(sig, ...
    sqrt((ols_s^2 * p^2 + sig^2 * n) / (p^2 + n)));
end % function nlrobustfit

%----------------------- Robust estimate of sigma
function s = madsigma(r,p)
%MADSIGMA    Compute sigma estimate using MAD of residuals from 0
n = length(r);
rs = sort(abs(r));
s = median(rs(max(1,min(n,p)):end)) / 0.6745;
end % function madsigma

% ---------------------- Jacobian
function J = getjacobian(beta,fdiffstep,model,X,yfit,nans,sweights)
    function yplus = call_model_nested(betaNew)
        yplus = model(betaNew, X);
        yplus(nans) = [];
    end

J = nirs.math.NonLinearModel.statjacobian(@call_model_nested, beta, fdiffstep, yfit(~nans));
% 
% [U,S,V]=svd(J,0);
% tol=sqrt(eps(norm(J,1)));
% J=U*(S+tol*eye(size(S)))*V';

if ~isempty(sweights)
    sweights = sweights(~nans);
    J = bsxfun(@times,sweights(:),J);
end
end % function getjacobian

%----------------------- Parse input arguments
function [errModel, weights, errModelParameter, options, iterative, maxweight] = parseInVarargin(varargin)
% Parse input arguments

argsToParse = varargin{:};
options = statset('nlinfit');
iterative = false;

% Process PVP - ErrorModel, Weights and ErrorModelParameter inputs should
% be passed in as PV pairs. statset can be passed in either as PV pair or
% directly as a struct.
numargs = numel(argsToParse);
if numargs > 0 && rem(numargs,2)
    % statset supplied directly as a struct.
    if isstruct(argsToParse{1}) || isempty(argsToParse{1})
        options = argsToParse{1};
        argsToParse = argsToParse(2:end);
    end
end
    
% Parse PV pairs
pnames = {'errormodel','weights', 'options', 'errorparameters', 'maxweight'};
defval = {'constant', 1, options, [], realmax};
[errModel, weights, options, errModelParameter, maxweight] = ...
    internal.stats.parseArgs(pnames, defval, argsToParse{:});

% Validate property values
if ~ischar(errModel)
    error(message('stats:nlinfit:InvalidErrorModel'));
end
ok = {'constant', 'proportional', 'exponential', 'combined'};
okv = find(strncmpi(errModel, ok, numel(errModel)));
if numel(okv) ~= 1
    error(message('stats:nlinfit:InvalidErrorModel'));
end
errModel = ok{okv};

if ~isnumeric(weights) && ~isa(weights, 'function_handle')
    error(message('stats:nlinfit:InvalidWeights'));
end

options = statset(statset('nlinfit'),options);

if numel(errModelParameter)>2 || ~isnumeric(errModelParameter)
    error(message('stats:nlinfit:BadErrorParam'))
end
switch errModel
    case 'combined'
        if isempty(errModelParameter)
            errModelParameter = [1 1];
        elseif numel(errModelParameter)~= 2
            % For combined error model, ErrorModelParameter should be a vector [a b]
            error(message('stats:nlinfit:BadCombinedParam', errModel));
        end
    case 'proportional'
        % Only a should be specified.
        if isempty(errModelParameter)
            errModelParameter = 1;
        elseif numel(errModelParameter)~=1
            error(message('stats:nlinfit:BadErrorParam1', errModel))
        end
    case {'constant', 'exponential'}
        % Only b should be specified.
        if isempty(errModelParameter)
            errModelParameter = 1;
        elseif numel(errModelParameter)~=1
            error(message('stats:nlinfit:BadErrorParam1', errModel))
        end
end

if ~isscalar(maxweight) || ~isreal(maxweight) || maxweight<=0
    error(message('stats:nlinfit:InvalidMaxWeight'));
end

% Check for conflicting error model and weights
if ~strcmpi(errModel, 'constant')
    if isa(weights, 'function_handle') || ~isscalar(weights) || weights~=1
        error(message('stats:nlinfit:ErrorModelWeightConflict'));
    end
end

% Robust fitting and weights
if ~isempty(options.RobustWgtFun) && (~strcmpi(errModel, 'constant') || isa(weights, 'function_handle') || ~isscalar(weights) || weights~=1)
    error(message('stats:nlinfit:ErrorModelRobustConflict'));
end

if any(strcmpi(errModel, {'proportional', 'combined'})) || isa(weights, 'function_handle')
    % Iteratively reweighted fitting required for proportional and
    % combined error model and weights that are a function of
    % predicted values
    iterative = true;
end

end % function parseInVarargin

%----------------------- Check weights
function nanweights = checkWeights(weights, yfit, y)
nanweights = zeros(size(y));

if (isnumeric(weights) && (~isscalar(weights) || weights~=1)) || isa(weights, 'function_handle')
    % If weights are set
    
    if isa(weights, 'function_handle')
        % function handle
        try
            wVec = weights(yfit);
        catch ME
            if isa(weights, 'inline')
                m = message('stats:nlinfit:InlineWeightFunctionError');
                throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
            elseif strcmp('MATLAB:UndefinedFunction', ME.identifier) ...
                    && ~isempty(strfind(ME.message, func2str(weights)))
                error(message('stats:nlinfit:WeightFunctionNotFound',func2str(weights)));
            else
                m = message('stats:nlinfit:WeightFunctionError',func2str(weights));
                throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
            end
        end
    else
        wVec = weights; % fixed weights
    end
    
    nanweights = isnan(wVec);
    % w should be real positive vector of the same size as y
    if any(~isreal(wVec)) || any(wVec(~nanweights)<=0) || numel(wVec)~=numel(y) ||~isvector(wVec) || ~isequal(size(wVec), size(y))
        error(message('stats:nlinfit:InvalidWeights'));
    end
    
end

end % function validateWeights

function e = error_ab(ab,y,f)
g = abs(ab(1)) + abs(ab(2))*abs(f);
e = sum(0.5*((y-f)./g).^2 + log(g));
end % function error_ab

function [y, model] = applyLogTransformation(y, model)
% Exponential, y = f*exp(a*e), or log(y) = log(f) + a*e

if ~isempty(y)
    if ~all(y>0)
        error(message('stats:nlinfit:PositiveYRequired'));
    else
        y = log(max(y,realmin));
    end
end

if ~isempty(model)
    % Exponential error model. Linearize the model as
    %   y = f*exp(a*e), or log(y) = log(f) + a*e
    model = @(phi,X) log(max(model(phi,X),realmin));
end

end % function applyModelTransformations

