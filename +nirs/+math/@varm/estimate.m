function [Mdl,standardErrors,logL,residuals,varargout] = estimate(Mdl,Y,varargin)
%ESTIMATE Estimate vector autoregression (VAR) models
%
% Syntax:
%
%	[EstMdl,EstSE,logL,E] = estimate(Mdl,Y)
%	[EstMdl,EstSE,logL,E] = estimate(Mdl,Y,name1,value1,...)
%
% Description:
%
%   Given an observed multivariate response series, estimate the parameters
%   of a VAR(p) model by maximum likelihood.
%
% Input Argument:
%
%   Mdl - VAR(p) model to fit to the multivariate time series Y (see below)
%     created by the VARM function.
%
%   Y - Multivariate response series to which the VAR(p) model Mdl is fitted. 
%     Y is a numobs-by-numseries matrix representing numobs observations of 
%     a numseries-dimensional time series. The last observation is the most 
%     recent.
%
% Optional Input Parameter Name/Value Pairs:
%
% 'Y0'      Presample responses to initialize the model. Y0 must be a matrix 
%           with the same number of columns (numseries) as Y. The number of 
%           observations (rows) of Y0 must be at least Mdl.P. If Y0 has more 
%           observations than necessary, only the most recent observations
%           are used. The last observation is the most recent. The default 
%           is to strip the first Mdl.P observations from the beginning of 
%           the responses Y, reducing the effective sample size.
%
% 'X'       Predictor data corresponding to a regression component in the 
%           model. X is a matrix with numpreds columns representing a
%           numpreds-dimensional time series of predictors. The number of
%           observations (rows) of X required depends upon the specification 
%           of presample responses Y0 (see above). If presample responses Y0 
%           are also specified, then the number of observations of X must be 
%           at least numobs; otherwise the number of observations of X must 
%           be at least numobs - Mdl.P to account for presample stripping. 
%           If X has more observations than necessary, only the most recent 
%           observations are used. The last observation is the most recent. 
%           The default is an empty series and the estimated VAR model has 
%           no regression component.
%
%'Display'  String or character vector indicating what estimation information 
%           to display in the command window. Choices are 'off' to display 
%           no information, 'table' to print a table of parameter estimates, 
%           standard errors, t statistics, and p-values, and 'full' to print 
%           all available estimation results. The default is 'off'.
%
%'MaxIterations'  Positive, integer maximum number of solver iterations (see 
%                 MVREGRESS for details). The default is 1000.
%
% Output Arguments:
%
%   EstMdl - Estimated VAR(p) model. 
%
%   EstSE - Estimated model parameter asymptotic standard errors. EstSE is 
%     a structure with the following parameter fields that correspond directly
%     to those in the estimated model (EstMdl):
%
%        - Constant
%        - AR (corresponding to lags 1,2,...,P)
%        - Beta
%        - Trend
%
%     The standard errors associated with any parameters held fixed as 
%     equality constraints are zero.
%
%   logL - Optimized loglikelihood objective function value.
%
%	  E - Residuals from the model fit. If presample responses Y0 are specified
%     explicitly, then E is a numobs-by-numseries matrix the same size as Y; 
%     otherwise the number of observations (rows) in E is numobs - Mdl.P to 
%     account for presample stripping.
%
% Notes:
%
%   o Observations with missing (NaN) values in the responses Y and predictors
%     X are removed via listwise deletion. That is, responses and predictors
%     are combined into a single time series, and all observations with a
%     NaN are removed. Similarly, missing values in the presample responses 
%     Y0 are also removed by listwise deletion.
%
% References:
% 
%   [1] Hamilton, J. D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994.
%
%   [2] Johansen, S. Likelihood-Based Inference in Cointegrated Vector
%       Autoregressive Models. Oxford: Oxford University Press, 1995.
%
%   [3] Juselius, K. The Cointegrated VAR Model. Oxford: Oxford University 
%       Press, 2006.
%
%   [4] Lutkepohl, H. New Introduction to Multiple Time Series Analysis.
%       Springer-Verlag, 2007.
%
% See also VARM.

% Copyright 2019 The MathWorks, Inc.

%
% Undocumented Usage: Weighted Estimation & Initial Model Estimates
%
% o Weighted Estimation:
%
%   To perform weighted estimation, VARM/ESTIMATE also accepts the following
%   name-value pair:
%
%  'Weights' Column vector of nonnegative weights to perform weighted estimation.
%            The length of Weights required depends upon the specification of 
%            presample responses Y0 (see above). If presample responses Y0 are
%            specified, then the length of Weights must be at least numobs + P;
%            otherwise the length of Weights must be at least numobs. If Weights
%            has more elements than necessary, only the most recent are used. 
%            The last weight is applied to the most recent observations of 
%            responses and predictors. The default is a vector of ones.
%
% o Initial Model Estimates:
%
%   To specify a VARM model that contains initial parameter estimates
%   subsequently passed to MVREGRESS, VARM/ESTIMATE also accepts the following
%   name-value pair:
%
%   'Mdl0'   Fully-specified VARM model containing initial parameter estimates.
%            Mdl0 must be of the same parametric form as the model to be 
%            estimated (Mdl). The parameter values in Mdl0 are passed to 
%            MVREGRESS as initial estimates (see optional 'beta0' and 'covar0'
%            inputs of MVREGRESS for details).
%

if numel(Mdl) > 1
   error(message('econ:varm:estimate:NonScalarModel'))
end

%
% Initialize some basic parameters and perform error-checking.
%

K = size(Mdl.PrivateConstant,1);     % # of response series (numseries) in Y(t)
P = Mdl.P;                           % # of pre-sample responses need for initialization

parser = inputParser;
parser.addRequired ('Y'            ,       @(x) validateattributes(x, {'double'}        , {'ncols' K 'nonempty'}, '', 'Responses'));
parser.addParameter('X'            ,   [], @(x) validateattributes(x, {'double'}        , {'2d'}                , '', 'Predictors'));
parser.addParameter('Y0'           ,   [], @(x) validateattributes(x, {'double'}        , {'ncols' K }          , '', 'Presample Responses'));
parser.addParameter('Display'      ,'off', @(x) validateattributes(x, {'char'  'string'}, {'scalartext'}        , '', 'Display'));
parser.addParameter('MaxIterations', 1000, @(x) validateattributes(x, {'double'}        , {'scalar' 'positive' 'integer'}, '', 'maximum number of iterations'));
parser.addParameter('Weights'      , ones(size(Y,1)+P,1), @(x) validateattributes(x, {'double'}, {'ncols' 1 'nonempty' 'nonnegative'}, '', 'Weights'));
parser.addParameter('Mdl0'         ,   [], @(x) validateattributes(x, {'nirs.math.varm'}          , {'scalar'}            , '', 'initial VARM model'));

parser.parse(Y, varargin{:});
   
Y             = parser.Results.Y;
X             = parser.Results.X;
Y0            = parser.Results.Y0;
display       = validatestring(lower(parser.Results.Display), {'off' 'table' 'full'}, '', 'Display');
maxIterations = parser.Results.MaxIterations;
W             = parser.Results.Weights;
Mdl0          = parser.Results.Mdl0;

%
% Set some flags and scrub the responses Y, predictors X, presample responses
% Y0, and weights W.
%
% The following code segment removes missing observations (NaNs) via listwise 
% deletion, and then ensures the presample responses Y0, predictors X, and
% weights W are of sufficient length.
%

updateDescription       =  Mdl.Description == getModelSummary(Mdl);                      % Should the description be updated?
isY0Specified           = ~any(strcmpi('Y0', parser.UsingDefaults));                     % Did the user specify presample responses Y0?
isRegressionIncluded    = ~any(strcmpi('X' , parser.UsingDefaults)) && (size(X,2) > 0);  % Did the user specify X with at least 1 predictors?
isInitialModelSpecified = ~any(strcmpi('Mdl0', parser.UsingDefaults));                   % Did the user specify initial estimated model?

%
% Note on Weighted Estimation:
%
% The following code segment will correctly handle missing data found in the
% weights W during the estimation period (i.e., the "in-sample" period common
% to Y, X, and W used for estimation). However, when optional presample Y0 
% is specified, it will NOT handle missing data in the presample period common
% to Y0 and W.
%

if any(isnan(X(:))) || any(isnan(Y(:))) || any(isnan(W(:)))
   [X,Y,W] = internal.econ.LagIndexableTimeSeries.listwiseDelete(X,Y,W);
end

if isY0Specified
   if any(isnan(Y0(:)))
      Y0 = internal.econ.LagIndexableTimeSeries.listwiseDelete(Y0);
   end
   try
      Y0 = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(P,K), 'Y0', Y0, P);
   catch exception
      exception.throwAsCaller();
   end
   T = size(Y,1);         % Sample size taken directly from responses Y
   W = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(size(Y,1) + P,1), 'Weights', W, size(Y,1) + P);
else
   T = size(Y,1) - P;     % Effective sample size to account for auto-stripping of presample responses
   W = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(size(Y,1),1), 'Weights', W, size(Y,1));
end

% 
% Ensure consistency between the predictor coefficients of the model and the
% number of predictors found in the input data X, and format the cell array
% of predictor information in block-diagonal form.
%

nPredictors = size(X,2);  % # of predictor series (X(t) data takes precedence)
I           = eye(K);

if isRegressionIncluded
   try
      X = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(T,size(X,2)), 'X', X, T);
   catch exception
      if strcmp(exception.identifier, 'econ:LagIndexableTimeSeries:checkPresampleData:InvalidPreSampleRows')
         error(message('econ:varm:estimate:InconsistentX',T))
      else
         exception.throwAsCaller();
      end
   end
   if isempty(Mdl.PrivateBeta)
      Mdl.PrivateBeta = nan(K,nPredictors);
   else
      if nPredictors ~= size(Mdl.PrivateBeta,2)
         error(message('econ:varm:estimate:InconsistentRegressionSpecification', nPredictors, size(Mdl.PrivateBeta,2)))
      end
   end
%
%  Format the cell vector of predictors in block-diagonal form consistent 
%  with the approach used below by MVREGRESS.
%
   XX = cell(T,1); 
   for t = 1:T
       XX{t} = kron(X(t,:), I);
   end

end

if isY0Specified
   Y = [Y0 ; Y];       % Prepend presample responses for convenience
end

%
% Set up parameter and parameter mapping vectors, and store the indices of
% each parameter included in the model. 
%
% Also, determine the number of model parameters included in the conditional 
% mean, which excludes parameters of the innovations covariance matrix which 
% are handled separately by MVREGRESS. The total # of parameters in the model 
% includes:
%
%    K constants
%    K^2 coefficients for each non-zero AR coefficient at positive lags (not necessarily = P)
%    K coefficients for each predictor
%    K trends
%

Lags        = Mdl.PrivateARLagOp.Lags;          % Lags included in the AR polynomial (includes zero)
nAR         = sum(Lags > 0);                    % # of AR lags included (not necessarily = P)
Lags        = Lags(Lags > 0);                   % Lags included in the difference equation (excludes zero)
nParameters = K + nAR*K*K + nPredictors*K + K;  % Total # of parameters included in the model (excludes innovations covariance matrix)
beta        = nan(nParameters,1);	              % Parameter vector (an element for each parameter included in the model)
solve       = true(nParameters,1);              % Parameter solve vector (TRUE = estimate, FALSE = do NOT estimate)
beta0       = zeros(nParameters,1);             % Initial parameter  estimates for MVREGRESS
covariance0 = cov(Y);                           % Initial covariance estimates for MVREGRESS

iConstant        = 1:K;
beta(iConstant)  = Mdl.PrivateConstant;
solve(iConstant) = isnan(beta(iConstant));
iStart           = iConstant(end) + 1;          % Starting index of next element

if isInitialModelSpecified
   beta0(iConstant) = Mdl0.PrivateConstant;
   covariance0      = Mdl0.PrivateCovariance;
end

if P > 0               % Are there any auto-regressive lags included?
   iFirstAR = iStart;  % Save the first index associated with the AR coefficients
   for i = Lags
       indices        =  iStart:(iStart - 1) + K*K;
     		beta(indices)  = -Mdl.PrivateARLagOp.Coefficients{i}(:);  % Reflect the coefficient
      	solve(indices) =  isnan(beta(indices));
       iStart         =  indices(end) + 1;
       if isInitialModelSpecified
          beta0(indices) = -Mdl0.PrivateARLagOp.Coefficients{i}(:);  % Reflect the coefficient
       end
   end
end

if isRegressionIncluded
   iBeta        = iStart:(iStart - 1) + K*nPredictors;
   beta(iBeta)  = Mdl.PrivateBeta(:);
   solve(iBeta) = isnan(beta(iBeta));
   iStart       = iBeta(end) + 1;
   if isInitialModelSpecified
      beta0(iBeta) = Mdl0.PrivateBeta;
   end
else
   iBeta = []; % Ensure it exists
end

iTrend        = iStart:(iStart - 1) + K;
beta(iTrend)  = Mdl.PrivateTrend;
solve(iTrend) = isnan(Mdl.PrivateTrend);
if isInitialModelSpecified
   beta0(iTrend) = Mdl0.PrivateTrend;
end

%
% Ensure that the innovations covariance matrix is not partially specified:
%
%  o It must either be all NaNs (completely unspecified), or 
%  o It must have no NaNs (completely specified).
%

if (sum(isnan(Mdl.PrivateCovariance(:))) > 0) && (sum(isnan(Mdl.PrivateCovariance(:))) < (K * K))
   error(message('econ:varm:estimate:InvalidCovariance'))
end

%
% Set a flag to determine if the covariance matrix is specified by the user. 
%

isCovarianceUnspecified = sum(isnan(Mdl.PrivateCovariance(:))) == (K * K);  % All NaNs?

%
% If the model is fully-specified, then simply compute the residuals and
% log-likelihood, set the fitted information to indicate that no estimation
% is performed, and return.
%

if all(~solve)
%
%  In the event the user has set all model parameters EXCEPT the innovations 
%  covariance matrix, calling INFER will throw an error. To prevent the error, 
%  the following code temporarily overwrites the covariance matrix, infers 
%  residuals, computes the covariance matrix from residuals, then infers 
%  the log-likelihood.
%
   if isCovarianceUnspecified
      Mdl.PrivateCovariance = eye(K);       % Overwrite if unspecified (avoids INFER error)
   end
   if isRegressionIncluded
      residuals = infer(Mdl,Y,'X',X);
      if isCovarianceUnspecified
         Mdl.PrivateCovariance = (residuals' * residuals)/T;
      end
      [~,logL]  = infer(Mdl,Y,'X',X);
   else
      Mdl.PrivateBeta = nan(K,0);           % Ensure model is consistent with X(t) data
      residuals       = infer(Mdl,Y);
      if isCovarianceUnspecified
         Mdl.PrivateCovariance = (residuals' * residuals)/T;
      end
      [~,logL] = infer(Mdl,Y);
      if updateDescription                  % Update model description if not user-specified
         Mdl.Description = getModelSummary(Mdl);
      end
   end
   standardErrors = getStandardErrors(Mdl,zeros(nParameters));
                                     
   Mdl.FitInformation.IsEstimated.Constant   = false(K,1);
   Mdl.FitInformation.IsEstimated.AR         = repmat({false(K)},1,P);
   Mdl.FitInformation.IsEstimated.Trend      = false(K,1);
   Mdl.FitInformation.IsEstimated.Beta       = false(K,nPredictors);
   Mdl.FitInformation.SampleSize             = T;
   Mdl.FitInformation.StandardErrors         = standardErrors;
   Mdl.FitInformation.LogLikelihood          = logL;
   Mdl.FitInformation.NumEstimatedParameters = 0;
   Mdl.FitInformation.Sigma = struct('Constant', zeros(K)                                , ...
                                     'AR'      , zeros(nParameters - K*(2 + nPredictors)), ...
                                     'Beta'    , zeros(K*nPredictors)                    , ...
                                     'Trend'   , zeros(K));
   if nargout > 4
      varargout = {Mdl.FitInformation.Sigma};
   end
   return
end

%
% Set up the multivariate regression problem.
%

D = cell(T,1);                            % Design cell array
R = Y(P+1:end,:);                         % In-sample responses stripped of pre-sample values
Z = zeros(K,nParameters);

Z(:,iConstant) = I;                       % Pull out of the FOR loop since it's deterministic

for t = (P + 1):(P + T)

    if P > 0               % Are there any auto-regressive lags included?
       iStart = iFirstAR;  % Restore the first index associated with the AR coefficients
   	   for i = Lags
           iColumns      = iStart:(iStart - 1) + K*K;
           Z(:,iColumns) = kron(Y(t-i,:),I);
           iStart        = iColumns(end) + 1;
       end
    end

   	if isRegressionIncluded
		     Z(:,iBeta) = XX{t-P};
    end

    Z(:,iTrend) = (t-P) * I;              % Weight the trend
%
%   Update the design array and adjust the responses to impose parameter
%   equality constraints.
%
    D{t-P}   =  sparse(Z(:,solve).*W(t));
    R(t-P,:) = sparse((R(t-P,:) - (Z(:,~solve)*beta(~solve))').*W(t));

end

%
% Estimate the parameters in the model.
%

sigma = zeros(nParameters);

[beta(solve), covariance, residuals, ...
 sigma(solve,solve), logL] = nirs.math.mvregress(D, R, 'covtype' , 'full'      , 'varformat', 'beta'       , ...
                                             'vartype' , 'fisher'    , 'maxiter'  , maxIterations, ...
                                             'beta0'   , beta0(solve), 'algorithm', 'ecm'        , ...
                                             'covar0',   covariance0);

% Record the pattern of NaNs in each model parameter (NaNs indicate estimated 
% parameters and non-NaNs indicate subset-constrained parameters held fixed).
%

Mdl.FitInformation.IsEstimated.Constant = isnan(Mdl.PrivateConstant);
Mdl.FitInformation.IsEstimated.AR       = mat2cell(isnan([Mdl.AR{:}]), (P > 0)*K, K(ones(1,P)));
Mdl.FitInformation.IsEstimated.Trend    = isnan(Mdl.PrivateTrend);
Mdl.FitInformation.IsEstimated.Beta     = isnan(Mdl.PrivateBeta);

%
% Unpack the parameter estimates from the vector returned by MVREGRESS and 
% pack them into the output VAR model.
%

Mdl.PrivateConstant = beta(iConstant);

if P > 0                       % Are there any auto-regressive lags included?
   iStart = iFirstAR;          % Restore the first index associated with the AR coefficients
	  for i = Lags
       indices                            =  iStart:(iStart - 1) + K*K;
       Mdl.PrivateARLagOp.Coefficients{i} = -reshape(beta(indices),K,K);    % Reflect the coefficient
       iStart                             =  indices(end) + 1;
   end
   iAR = iFirstAR:indices(end);
else
   iAR = [];  % Ensure it exists
end

if isRegressionIncluded
   Mdl.PrivateBeta = reshape(beta(iBeta), K, nPredictors);
else
%
%  The presence of X(t) takes precedence, so overwrite any regression
%  coefficients to ensure the model is consistent.
%
   Mdl.PrivateBeta = zeros(K,0);
end

Mdl.PrivateTrend = beta(iTrend);

%
% If the innovations covariance matrix is unspecified, then it acts as a
% derived parameter and by default we simply assign the sample residual 
% covariance obtained from MVREGRESS and report the corresponding log-likelihood.
%
% However, if the covariance matrix is specified by the user, then this
% introduces a discrepancy between it and the sample residual covariance and
% so we must re-evaluate the log-likelihood.
%

if isCovarianceUnspecified
   Mdl.PrivateCovariance = covariance;   % Innovations covariance from MVREGRESS
else
   if isRegressionIncluded
      [~,logL] = infer(Mdl,Y,'X',X);
   else
      Mdl.PrivateBeta = nan(K,0);        % Ensure model is consistent with X(t) data
      [~,logL] = infer(Mdl,Y);
   end    
end

%
% Format the structure of parameter standard errors.
%

standardErrors = getStandardErrors(Mdl,sigma);

%
% Pack estimation results into the model for summary display.
%

Mdl.FitInformation.SampleSize             = T;
Mdl.FitInformation.StandardErrors         = standardErrors;
Mdl.FitInformation.LogLikelihood          = logL;
Mdl.FitInformation.NumEstimatedParameters = sum(solve);
Mdl.FitInformation.Y0                     = Y(1:P,:);
Mdl.FitInformation.Sigma                  = struct('Constant', sigma(iConstant,iConstant), ...
                                                   'AR'      , sigma(iAR,iAR)            , ...
                                                   'Beta'    , sigma(iBeta,iBeta)        , ...
                                                   'Trend'   , sigma(iTrend,iTrend));
%
% Calling this method with more than 4 outputs is an indication that the
% error covariance matrices associated with each parameter of the VARX model  
% are needed. Specifically, it indicates that subsequent calculations involve 
% the computation of asymptotic standard errors of parameters associated 
% with the Johansen (i.e., the first) step of VEC model estimation.
%
% Notice that the entire error covariance matrix is NOT returned, but rather 
% a data structure whose fields are packed with the error covariance "blocks" 
% of each parameter in the model. This makes subsequent manipulation much more
% convenient.
%
% This is an UNDOCUMENTED feature, and is only needed when this method is
% called by the VECM/ESTIMATE method for cointegration rank 1 <= r <= (K-1).
%
% External users should NOT rely on this behavior!
%

if nargout > 4
   varargout = {Mdl.FitInformation.Sigma};
end

%
% Update the description of the estimated model ONLY if the description of 
% the input model was auto-generated. Otherwise, retain a user-specified 
% description. 
%
% The description is deemed to be auto-generated if the description in the 
% original input model matches the model summary string.
%

if updateDescription
   Mdl.Description = getModelSummary(Mdl);
end

%
% Display summary information if requested.
%

if strcmpi(display, 'full')
   summarize(Mdl);
elseif strcmpi(display, 'table')
   tbl = summarize(Mdl); 
   disp(tbl.Table)
end

end   % End of ESTIMATE method

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function standardErrors = getStandardErrors(Mdl,sigma)
%
% Map the row/column indices of the parameter covariances of the estimated
% VARX model to the standard errors of individual coefficients. The results
% are packed into a data structure with the same fieldnames as the corresponding 
% properties of the VAR model.
%

K           = size(Mdl.PrivateConstant,1);  % # of response series (numseries) in Y(t)
P           = Mdl.P;                        % # of pre-sample responses need for initialization
nPredictors = size(Mdl.PrivateBeta,2);      % # of predictors in X(t)

Lags        = Mdl.PrivateARLagOp.Lags;      % Lags included in the AR polynomial (includes zero)
nAR         = sum(Lags > 0);                % # of AR lags included (not necessarily = P)
Lags        = Lags(Lags > 0);               % Lags included in the difference equation (excludes zero)

iConstant   = 1:K;                          % Indices of VARX model constant
iNext       = iConstant(end) + 1;           % Starting index of next element

iAR = iNext:(iNext - 1 + nAR * K * K);
if nAR > 0
   iNext = iAR(end) + 1;
end

iBeta = iNext:(iNext - 1 + K * nPredictors);
if nPredictors > 0
   iNext = iBeta(end) + 1;
end

iTrend = iNext:(iNext - 1 + K);

%
% Compute the standard errors of all reported model parameters.
%

standardErrors.Constant = reshape(sqrt(diag(sigma(iConstant,iConstant))),K,1);

if P > 0
   iAR = reshape(iAR,K*K,nAR);
   for i = 1:P
       if any(i == Lags)
          iLag = (i == Lags);
          standardErrors.AR{i} = reshape(sqrt(diag(sigma(iAR(:,iLag),iAR(:,iLag)))),K,K);
       else
          standardErrors.AR{i} = zeros(K);
       end
   end
else
   standardErrors.AR = {};
end

standardErrors.Beta  = reshape(sqrt(diag(sigma(iBeta ,iBeta))) ,K,nPredictors);
standardErrors.Trend = reshape(sqrt(diag(sigma(iTrend,iTrend))),K,1);

%
% To assign the full error covariance matrix (advanced users only!), uncomment
% the following line of code.
%
% standardErrors.EstParamCov = sigma;

end
