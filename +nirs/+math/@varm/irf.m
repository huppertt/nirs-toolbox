function [baseline,lowerBounds,upperBounds,IRF] = irf(Mdl,varargin)
%IRF Impulse response function of vector autoregression (VAR) models
%
% Syntax:
%
%   [response,lower,upper] = irf(Mdl)
%   [response,lower,upper] = irf(Mdl,name1,value1,name2,value2,...)
%
% Description:
%
%   Compute the dynamic response, or the impulse response function (IRF), to
%   a one standard deviation innovation shock to each variable in the system. 
%   Optionally, the lower and upper confidence bounds of the IRF are also 
%   computed by Monte Carlo simulation or bootstrapping. 
%
% Input Arguments:
%
%   Mdl - VAR model created by the VARM constructor or VARM/ESTIMATE method.
%
% Optional Input Parameter Name/Value Pairs:
%
% 'NumObs'  Positive integer indicating the number of observations included
%        in the IRF (the number of periods for which the impulse response is 
%        computed). The default is 20.
%
% 'Method'  String or character vector indicating the method by which the 
%        IRF is computed. Options are 'orthogonalized' and 'generalized'. 
%        The default is 'orthogonalized'.
%
% Optional Input Parameter Name/Value Pairs Used Only for Confidence Bounds:
%
% 'NumPaths'  Positive integer indicating the number of sample paths (trials).
%        The default is 100.
%
% 'SampleSize'  Positive integer indicating the number of observations 
%        simulated or bootstrapped per path. If the model Mdl is estimated, 
%        the default number of observations is the sample size used for 
%        estimation (see VARM/SUMMARIZE). If Mdl is not estimated and 
%        confidence bounds are bootstrapped, the default sample size is the 
%        number of observations (rows) in the residual series (E) used for 
%        bootstrapping. If the model Mdl is not estimated and confidence 
%        bounds are simulated by Monte Carlo, the sample size must be 
%        specified. 
%
% 'Y0'   Presample responses to initialize the model. Y0 must be a matrix 
%        with the same number of columns (numseries) as the number of series
%        found in Mdl (Mdl.NumSeries). The number of observations (rows) of 
%        Y0 must be at least Mdl.P. If Y0 has more observations than necessary,
%        only the most recent observations are used. The last observation is
%        the most recent. If the model Mdl is estimated, the default Y0 are 
%        those used for estimation. If the model Mdl is not estimated and 
%        confidence bounds are computed, Y0 is required. 
%
% 'X'    Predictor data corresponding to a regression component in the model. 
%        X is a 2-D matrix with numpreds columns (a numpreds-dimensional time 
%        series of predictors). The last observation is the most recent. The 
%        number of observations in X must equal or exceed the sample size (see
%        'SampleSize'). When the number of observations in X exceeds the number
%        necessary, only the most recent observations are used. By default 
%        the model has no regression component. 
%
% 'E'    A numperiods-by-numseries matrix of residuals from which to bootstrap
%        (resample with replacement) confidence bounds. The default computes 
%        confidence bounds by Monte Carlo simulation. 
%
% 'Confidence'  Nonnegative scalar confidence level. Confidence (C) is such 
%        that 100*(1-C)/2 percent of the impulse responses lie below and above 
%        the lower and upper confidence bounds, respectively. Confidence must 
%        be between 0 and 1 (inclusive). The default is 0.95 (95% confidence 
%        bounds).
%
% Output Arguments:
%
%   response - numobs-by-numseries-by-numseries array of impulse responses.
%       The IRF in element (t+1,i,j) is the response of variable j to a one 
%       standard deviation innovation shock to variable i at time 0 for 
%       t = 0, 1, ..., numobs-1.
%
%   lower - numobs-by-numseries-by-numseries array of lower confidence 
%       bounds of impulse responses.
%
%   upper - numobs-by-numseries-by-numseries array of upper confidence 
%       bounds of impulse responses.
%
% Notes:
%
%   o A model Mdl is estimated if it is obtained directly from the ESTIMATE 
%     function and remains unchanged thereafter. 
%
%   o The orthogonalized IRF method orthogonalizes the innovation shocks via
%     the Cholesky factorization of the model covariance matrix (Mdl.Covariance).
%     When orthogonalized, the impulse responses generally depend on the 
%     ordering of the variables. In contrast, the generalized IRF method is 
%     invariant to the ordering of the variables, and therefore the two methods
%     generally produce different results.
%
%     However, if Mdl.Covariance is diagonal, then the results of both methods
%     are identical. If Mdl.Covariance is non-diagonal, then the results will 
%     coincide only when the first variable is shocked (i.e., the first column
%     of each page will always coincide).
%
%   o Missing values, indicated by NaNs, are removed from Y0, X, and E by 
%     listwise deletion (any row of a series containing at least one NaN 
%     is removed), reducing the effective sample size. 
%
%   o The inclusion of a regression component is based on the presence of 
%     the predictor matrix X, which represents a single path of a multivariate
%     time series. When the model includes a regression component, the entire
%     predictor matrix X is applied to all paths used to compute confidence 
%     bounds.
%
% Example:
%
%   Fit a VAR(2) model to the Johansen Danish data and estimate 95% confidence 
%   bounds by bootstrapping 100 paths of the fitted residuals. Plot the IRF
%   and 95% confidence bounds of the bond rate (variable 3) to a shock to 
%   real income (variable 2).
%
%   load Data_JDanish
%   rng(1)
%   Y = Data;
%   [Fit,~,~,E] = estimate(varm(4,2), Y);
%   [response,lower,upper] = irf(Fit, 'E', E);
%   plot(0:19,[response(:,2,3) lower(:,2,3) upper(:,2,3)])
%   xlabel('Observation Time');
%   ylabel('Response');
%   grid
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
% See also ARMAIRF, SIMULATE, FILTER.

% Copyright 2018 The MathWorks, Inc.

if numel(Mdl) > 1
   error(message('econ:varm:irf:NonScalarModel'))
end

%
% Ensure the input model is fully-specified.
%

if any(isnan(Mdl.PrivateConstant))
   error(message('econ:varm:irf:UnspecifiedConstant'))
end

if any(any(isnan([Mdl.AR{:}])))
   error(message('econ:varm:irf:UnspecifiedAR'))
end

if any(isnan(Mdl.PrivateTrend))
   error(message('econ:varm:irf:UnspecifiedTrend'))
end

if any(any(isnan(Mdl.PrivateCovariance)))
   error(message('econ:varm:irf:UnspecifiedCovariance'))
end

%
% Initialize some basic parameters and perform error-checking.
%

K = size(Mdl.PrivateConstant,1);     % # of response series (numseries) in Y(t)
P = Mdl.P;                           % # of pre-sample responses need for initialization

parser = inputParser;
parser.addParameter('NumObs'    ,  20             , @(x) validateattributes(x, {'double'}       , {'scalar' 'integer' '>' 0}, '', 'number of observations'));
parser.addParameter('NumPaths'  ,  100            , @(x) validateattributes(x, {'double'}       , {'scalar' 'integer' '>' 0}, '', 'number of paths'));
parser.addParameter('SampleSize',  0              , @(x) validateattributes(x, {'double'}       , {'scalar' 'integer' '>' 0}, '', 'sample size'));
parser.addParameter('Method'    , 'Orthogonalized', @(x) validateattributes(x, {'char' 'string'}, {'scalartext'}            , '', 'impulse response method'));
parser.addParameter('Y0'        ,  []             , @(x) validateattributes(x, {'double'}       , {'2d' 'ncols'  K}         , '', 'presample responses'));
parser.addParameter('X'         ,  []             , @(x) validateattributes(x, {'double'}       , {'2d'}                    , '', 'predictors'));
parser.addParameter('E'         ,  []             , @(x) validateattributes(x, {'double'}       , {'2d' 'ncols'  K}         , '', 'bootstrapped residuals'));
parser.addParameter('Confidence',  0.95           , @(x) validateattributes(x, {'double'}       , {'scalar' '>=' 0 '<=' 1}  , '', 'confidence interval'));

parser.parse(varargin{:});

numObs     = parser.Results.NumObs;
nPaths     = parser.Results.NumPaths;
T          = parser.Results.SampleSize;
Y0         = parser.Results.Y0;
X          = parser.Results.X;
E          = parser.Results.E;
method     = validatestring(parser.Results.Method, {'orthogonalized', 'generalized'}, '', 'impulse response method');
confidence = parser.Results.Confidence;

%
% Set some flags for later use.
%

isY0Specified        = ~any(strcmpi('Y0', parser.UsingDefaults));   % Did the user specify presample responses Y0?
isRegressionIncluded = ~any(strcmpi('X' , parser.UsingDefaults));   % Did the user specify predictors X?
nPredictors          =  size(X,2);                                  % # of predictor series (X(t) data takes precedence)
isBootstrapping      = ~any(strcmpi('E' , parser.UsingDefaults));   % Are confidence bounds computed by bootstrapping?
isModelEstimated     = ~isempty(Mdl.FitInformation);                % Is the model estimated by the ESTIMATE method?

%
% Perform some basic error checking.
%

if (nargout > 1) && ~isModelEstimated && ~isY0Specified
   error(message('econ:varm:irf:UnspecifiedY0'))
end

if (nargout > 1) && ~isModelEstimated && ~isBootstrapping && any(strcmpi('SampleSize', parser.UsingDefaults))
   error(message('econ:varm:irf:UnspecifiedSampleSize'))
end

%
% Listwise-delete any missing observations from the presample responses Y0(t), 
% residuals E(t), and predictors X(t). 
%

if any(isnan(X(:)))    % Checking X(t)
   X = internal.econ.LagIndexableTimeSeries.listwiseDelete(X);
end

if any(isnan(E(:)))    % Checking E(t)
   E = internal.econ.LagIndexableTimeSeries.listwiseDelete(E);
end

if any(isnan(Y0(:)))   % Checking Y0(t)
   Y0 = internal.econ.LagIndexableTimeSeries.listwiseDelete(Y0);
end

%
% Determine the sample size (T) if unspecified.
%

if (nargout > 1) && any(strcmpi('SampleSize', parser.UsingDefaults))
   if isModelEstimated
      results = summarize(Mdl);
      T       = results.SampleSize;
   else
      T = size(E,1);
   end
end

%
% Determine the pre-sample responses (Y0) if unspecified.
%

if (nargout > 1) && ~isY0Specified && isModelEstimated
   Y0 = Mdl.FitInformation.Y0;
end

%
% Missing observations (NaNs) have been removed by listwise deletion, so now 
% ensure that they are of sufficient length. The following code segment also 
% ensures consistency between the predictor coefficients of the model and the 
% number of predictors in X(t).
%

if nargout > 1         % Do we compute confidence bounds?
   try
      Y0 = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(P,K), 'Y0', Y0, P);
   catch exception
      exception.throwAsCaller();
   end

   if isRegressionIncluded
      try 
         X = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(T,size(X,2)), 'X', X, T);
      catch exception
         if strcmp(exception.identifier, 'econ:LagIndexableTimeSeries:checkPresampleData:InvalidPreSampleRows')
            error(message('econ:varm:irf:InconsistentX',T))
         else
            exception.throwAsCaller();
         end
      end    
   
      if isempty(Mdl.PrivateBeta) || any(any(isnan(Mdl.PrivateBeta)))
         error(message('econ:varm:irf:UnspecifiedBeta'))
      else
         if nPredictors ~= size(Mdl.PrivateBeta,2)
            error(message('econ:varm:irf:InconsistentRegressionSpecification', nPredictors, size(Mdl.PrivateBeta,2)))
         end
      end
   end
end

%
% Compute the baseline IRF analytically using the utility ARMAIRF.
%

baseline = armairf(Mdl.AR, {}, 'InnovCov', Mdl.PrivateCovariance, 'NumObs', numObs, 'Method', method);

%
% Now compute the approximate lower/upper confidence bounds. In what follows, 
% we adopt the guidance outlined in Lutkepohl (Section 3.7.4 pages 126-129, 
% and Appendix D.3, pages 709-712, for bootstrapping).
%

if nargout > 1                      % Do we compute confidence bounds?

   toFit             = varm(K,P);   % Create a template without constraints.
   toFit.Description = "";          % Avoid updating description for performance.
%
%  Since VAR models do not, by default, include time trends and regression
%  components, pre-allocate those now and update later if necessary.
%
   if any(Mdl.PrivateTrend)         % Add an unconstrained linear time trend in case Mdl is NOT estimated.
      toFit.Trend = nan(K,1);       % If estimated, this is updated below.
   end
   toFit.Beta = nan(K,nPredictors); % If estimated, this is updated below.
%
%  Impose subset constraints if the original input model (Mdl) is estimated
%  and such information is available.
%
%  Any subset constraints included in the estimation of the original input 
%  model (Mdl) are reflected in the simulated/bootstrapped responses upon 
%  which subsequent estimations are based. Therefore, in many situations 
%  explicitly imposing the same subset constraints makes little difference
%  in the simulated/bootstrapped statistics when averaged over a large sample.  
%  However, in some situations the imposition of subset constraints does 
%  matter, and so we impose them below.
%
   if isModelEstimated
      toFit.PrivateConstant(~Mdl.FitInformation.IsEstimated.Constant) = Mdl.PrivateConstant(~Mdl.FitInformation.IsEstimated.Constant);
      toFit.PrivateTrend(~Mdl.FitInformation.IsEstimated.Trend)       = Mdl.PrivateTrend(~Mdl.FitInformation.IsEstimated.Trend);
      if isRegressionIncluded
         toFit.PrivateBeta(~Mdl.FitInformation.IsEstimated.Beta) = Mdl.PrivateBeta(~Mdl.FitInformation.IsEstimated.Beta);
      end
      for i = 1:P
          toFit.AR{i}(~Mdl.FitInformation.IsEstimated.AR{i}) = Mdl.AR{i}(~Mdl.FitInformation.IsEstimated.AR{i});
      end
   end

   if isBootstrapping
%
%     Bootstrap responses Y(t) as outlined in Lutkepohl Appendix D.3:
%       o Center the residuals
%       o Draw random observation times (row indices) 
%       o Resample residuals with replacement
%       o Filter resampled residuals to compute filtered responses
%
      E = E - mean(E);
      t = randi(size(E,1), [T nPaths]);
      Z = permute(reshape(E(t(:),:), [T nPaths K]), [1 3 2]);
      if isRegressionIncluded
         Y = filter(Mdl, Z, 'Y0', Y0, 'X', X, 'Scale', false);
      else
         Y = filter(Mdl, Z, 'Y0', Y0, 'Scale', false);
      end
   else
%
%     Simulate responses Y(t) via Monte Carlo.
%
      if isRegressionIncluded
         Y = simulate(Mdl, T, 'NumPaths', nPaths, 'Y0', Y0, 'X', X);
      else
         Y = simulate(Mdl, T, 'NumPaths', nPaths, 'Y0', Y0);
      end
   end
%
%  For each simulated/bootstrapped sample response path in Y(t), re-estimate 
%  the model using the same presample responses Y0(t) to ensure the same initial 
%  conditions. Then compute the IRF for each fitted model and accumulate IRF 
%  results.
%
%  In what follows, the path-by-path IRF results for all N paths are stored 
%  in the pages (3rd dimension) as follows:
%
%           Path #1                          Path #2                ...           Path #N
%  ------------------------------   ------------------------------  ...  ------------------------------  
%  IRF_1(t) IRF_2(t) ... IRF_K(t)   IRF_1(t) IRF_2(t) ... IRF_K(t)  ...  IRF_1(t) IRF_2(t) ... IRF_K(t)  
%
   IRF = nan(numObs, K, K * nPaths);
   
   for iPath = 1:nPaths
       if isRegressionIncluded
          Fit = estimate(toFit, Y(:,:,iPath), 'Y0', Y0, 'X', X);
       else
          Fit = estimate(toFit, Y(:,:,iPath), 'Y0', Y0);
       end
%
%      Compute and store the IRF for the current sample path.
%
       iPages          = ((iPath - 1)*K + 1):(iPath*K);
       IRF(:,:,iPages) = armairf(Fit.PrivateARLagOp, {}, 'InnovCov', Fit.PrivateCovariance, 'NumObs', numObs, 'Method', method);
   end
%
%  Compute lower/upper confidence bounds as quantiles of the IRF distribution.
%
   lowerBounds   = nan(size(baseline));
   upperBounds   = nan(size(baseline));
   probabilities = [(1 - confidence)/2  (1 + confidence)/2];

   for iResponse = 1:K         % The variable whose response to the shocked variable is measured
       for iShock = 1:K        % The shocked variable
           quantiles = quantile(IRF(:,iShock,iResponse:K:end), probabilities, 3);
           lowerBounds(:,iShock,iResponse) = quantiles(:,1);
           upperBounds(:,iShock,iResponse) = quantiles(:,2);
       end
   end

end

end