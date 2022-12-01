function [Y,YMSE] = forecast(Mdl, horizon, Y0, varargin)
%FORECAST Forecast vector autoregression (VAR) model responses
%
% Syntax:
%
%   [Y,YMSE] = forecast(Mdl,numperiods,Y0)
%   [Y,YMSE] = forecast(Mdl,numperiods,Y0,param1,val1,param2,val2,...)
%
% Description:
%
%   Forecast responses of a vector autoregression (VAR) model. 
%
% Input Arguments:
%
%   Mdl - VAR model created by the VARM constructor or VARM/ESTIMATE method.
%
%   numperiods - Positive integer specifying the number of periods in the 
%     forecast horizon. 
%
%   Y0 - Presample response matrix with numseries columns providing initial 
%     values for the forecast. When computing unconditional forecasts, Y0 is 
%     a numobs-by-numseries-by-numpaths 3-D matrix in which each page 
%     represents a single path of presample responses. When conditional 
%     forecasts are specified (see input YF below), Y0 is a 2-D or 3-D matrix. 
%     As a 2-D matrix, initial values in Y0 are applied to each path of Y 
%     and all output paths evolve from common initial conditions. As a 3-D 
%     matrix, Y0 provides initial conditions for each output path. For all 
%     pages observations across any row occur at the same time. Y0 may have 
%     any number of rows, provided at least Mdl.P observations exist to 
%     initialize the model. If the number of rows exceeds Mdl.P, then only 
%     the most recent Mdl.P observations are used. The last observation of 
%     any page is the most recent.
%
% Optional Input Parameter Name/Value Pairs:
%
%  'X'  Time series of forecasted (future) predictors to include a regression
%       component in the model forecast. X is a 2-D matrix with numpreds 
%       columns representing a numpreds-dimensional series of predictor 
%       forecasts, and the first row contains the 1-period-ahead forecast, 
%       the second row the 2-period-ahead forecast, and so forth. X may 
%       have any number of rows, provided at least numperiods rows exist 
%       to cover the forecast horizon. If the number of rows exceeds the 
%       forecast horizon, then only the first numperiods rows are used. By 
%       default the forecasts include no regression component.
%
%  'YF' Time series of forecasted (future) responses to support conditional 
%       forecasting. YF is a 2-D or 3-D matrix with numseries columns, and 
%       captures the state or knowledge of responses as they evolve from 
%       the presample past (Y0) into the future. The first row contains the
%       1-period-ahead forecast, the second row the 2-period-ahead forecast, 
%       and so forth. As a 2-D matrix, known values in YF are treated as 
%       deterministic forecasts and are applied to each of numpaths output 
%       paths. As a 3-D matrix, YF provides deterministic forecasts for each 
%       output path. YF may have any number of rows, provided at least 
%       numperiods rows exist to cover the forecast horizon. If the number 
%       of rows exceeds the forecast horizon, then only the first numperiods 
%       rows are used. By default, YF is an array of NaN values indicating 
%       a complete lack of knowledge of the future state of all responses, 
%       and the output responses Y (see below) are conventional minimum 
%       mean square error (MMSE) forecasts.
%
% Output Arguments:
%
%   Y - numperiods-by-numseries-by-numpaths matrix of minimum mean square 
%     error (MMSE) response forecasts. For all paths, the first row of Y 
%     contains the 1-period-ahead forecast, the second row contains the 
%     2-period-ahead forecast, and so forth. The last row contains the 
%     response forecast at the forecast horizon.
%
%   YMSE - numperiods-by-1 cell vector of mean square errors (MSE) of the
%     response forecasts in Y (i.e., a time series of numseries-by-numseries 
%     error covariance matrices). The first element of YMSE contains the 
%     1-period-ahead forecast error covariance, the second element contains 
%     the 2-period-ahead forecast error covariances, and so forth. The last
%     element contains the forecast error covariance at the forecast horizon.
%     YMSE is identical for all paths.
%
% Notes:
%
%   o When computing unconditional forecasts, the number of pages (numpaths)
%     in the output response forecasts Y is determined by the number of 
%     pages in the input presample responses Y0. However, when computing 
%     conditional forecasts, the number of pages in the output Y is determined
%     by the number of pages in the input presample responses Y0 and the input 
%     response forecasts YF. If both Y0 and YF inputs have more than one page, 
%     then the number of pages in the output Y is the smaller of the number 
%     of pages in Y0 and YF; otherwise, the number of pages in Y is the larger 
%     of the number of pages in the Y0 and YF inputs. If the number of pages
%     in the Y0 or YF inputs is greater than one and exceeds the number of 
%     output paths (numpaths), only the first numpaths pages are used.
%      
%   o Missing values, indicated by NaNs, are removed from Y0 by listwise 
%     deletion, reducing the effective sample size. That is, the 3-D presample
%     responses Y0 are first converted to an equivalent 2-D series by 
%     horizontal concatenation, and any row of the concatenated series with 
%     at least one NaN is then removed. This ensures that the forecasted 
%     responses of each path are the same size and based on the same forecast
%     times, yet also implies that in the presence of missing forecasts the 
%     results obtained from multiple paths of Y0 may differ from those 
%     obtained from each path individually.
%
%     Similarly, missing values in the predictor forecasts X are also removed
%     by listwise deletion. Additionally, when performing a conditional 
%     forecast (see below) the rows of predictor forecasts X and response 
%     forecasts YF are associated with the same future observation times, and 
%     any rows of YF associated with any NaNs in X are also removed.
%
%   o When YF is specified as an input, NaNs also indicate missing forecast 
%     values, but also indicate the pattern of partial knowledge of the state
%     of future responses. Non-missing values in YF serve as deterministic 
%     forecasts known in advance (e.g., set by policy), and in such cases any 
%     corresponding missing NaN values are forecasted conditional on 
%     non-missing information.
%
%     When missing information is included in the input response forecasts 
%     YF, the missing forecasts are derived from a Kalman filter. The input 
%     model is converted to its equivalent state-space representation in 
%     which the missing forecasts are the filtered states computed as a
%     conditional expectation (see SSM/FILTER for additional details).
%
%     When input YF is unspecified, or otherwise contains only missing values,
%     the output response forecasts Y(t) are conventional unconditional 
%     minimum mean square error (MMSE) forecasts.
%
%  o  The inclusion of a regression component is based on the presence of 
%     the predictor matrix X. Although each page of the output response
%     forecast Y represents a different path, the predictor matrix X 
%     represents a single path of a multivariate time series. When the 
%     model forecast includes a regression component, the entire predictor 
%     matrix X is applied to every path of the output response forecast Y.
%
%  o  Although only the first Mdl.P observations of presample responses Y0
%     are explicitly used to initialize the model, the total number of 
%     observations (with no missing values) determines the time origin 
%     associated with models that include linear time trends. 
%
%     The time origin (t0) used to compute the trend component is defined as 
%     the number of observations in Y0(t) adjusted for the number of lagged 
%     responses needed to initialize the model:
%
%     t0 = size(Y0,1) - Mdl.P
%
%     Therefore, the times used in the trend component of the output responses
%     Y(t) are based on times t = t0 + 1, t0 + 2, ..., t0 + numobs. This 
%     convention is consistent with the default behavior of model estimation
%     in which the first Mdl.P responses are stripped, reducing the effective
%     sample size.
%
%   o If the input model (Mdl) includes a regression component, then the
%     predictor forecast matrix X is typically specified to maintain model
%     consistency. For example, if the input model was estimated and included 
%     a regression component by specifying predictors X, then the forecasts 
%     of the same predictors are generally specified as well.
%
%   o When computing the mean square error of the response forecasts (YMSE), 
%     the predictor data (X) are treated as exogenous, non-stochastic, and
%     statistically independent of the model innovations. Therefore, the 
%     output YMSE reflects the error covariance associated with the 
%     autoregressive component of the input model Mdl alone.
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
% See also SIMULATE.

% Copyright 2017 The MathWorks, Inc.   

if numel(Mdl) > 1
   error(message('econ:varm:forecast:NonScalarModel'))
end

%
% Initialize some basic parameters and perform error-checking.
%

K = size(Mdl.PrivateConstant,1);   % # of responses (numseries) in Y(t)
P = Mdl.P;                         % # of pre-sample responses need for initialization

parser = inputParser;
parser.addRequired ('numperiods',     @(x) validateattributes(x, {'double'}, {'scalar' 'integer' '>' 0}, '', 'forecast horizon'));
parser.addRequired ('Y0'        ,     @(x) validateattributes(x, {'double'}, {'ncols'   K             }, '', 'presample responses'));
parser.addParameter('X'         , [], @(x) validateattributes(x, {'double'}, {'2d'}                    , '', 'forecasted predictor matrix'));
parser.addParameter('YF'        ,NaN, @(x) validateattributes(x, {'double'}, {'ncols'   K }            , '', 'forecasted responses'));
parser.parse(horizon, Y0, varargin{:});

Y0      = parser.Results.Y0;
X       = parser.Results.X;
Y       = parser.Results.YF;
horizon = parser.Results.numperiods;

%
% Determine whether an unconditional or conditional forecast is selected.
%
% A conditional forecast is indicated the presence of an input response 
% series Y(t) whose elements are NOT all NaNs.
%

if (any(strcmpi('YF', parser.UsingDefaults)) || all(isnan(Y(:))))
%
%  No response forecasts are specified, so compute the unconditional minimum 
%  mean square error (MMSE) forecast in the usual sense.
%
%  Remove optional parameters specific to FORECAST and call FILTER.
%
   varargin(1:2:end) =  cellstr(varargin(1:2:end));
   indices           = ~strcmpi('YF', varargin(1:2:end));
   varargin          =  reshape(varargin, 2, numel(varargin)/2);
   varargin          =  varargin(:,indices);
   nPaths            =  size(Y0,3);
%
%  Now that YF has been removed from the optional N-V pairs, the only optional
%  input that can remain in VARARGIN is the forecasted predictors X(t). 
%
%  Since FILTER will strip observations from the beginning of X(t), and FORECAST 
%  needs to strip the most distant observations, we obviate the distinction by 
%  retaining only the first "horizon" observations.
%
   if ~isempty(varargin)
      if any(isnan(varargin{end}(:)))    % Check X(t) for missing observations
         varargin{end} = internal.econ.LagIndexableTimeSeries.listwiseDelete(varargin{end});
      end
      if size(varargin{end},1) < horizon
         error(message('econ:varm:forecast:InconsistentX', horizon))
      else
         varargin{end} = varargin{end}(1:horizon,:);
      end
   end

   try
      Y = filter(Mdl, zeros(horizon,K,nPaths), 'Y0', Y0, varargin{:}, 'Scale', false);
   catch exception
      exception.throwAsCaller();
   end

%
%  Compute mean square error covariance if necessary.
%
   if nargout > 1
%
%     Compute forecast error covariances by converting the VAR model to its 
%     truncated infinite-degree MA representation.
%
      wState  = warning;                           % Save warning state
      cleanUp = onCleanup(@() warning(wState));    % Restore warning state
   
      warning('off', 'econ:LagOp:mldivide:WindowNotOpen')   
      warning('off', 'econ:LagOp:mldivide:WindowIncomplete')

      MA      = mldivide(Mdl.PrivateARLagOp, {eye(K)}, 'Degree', horizon - 1, 'RelTol', 0, 'AbsTol', 0);
      YMSE    = cell(horizon,1);
      YMSE{1} = Mdl.PrivateCovariance;   % Since MA{1} = I = eye(K)
   
      for t = 2:horizon
          YMSE{t} = YMSE{t-1} + MA.Coefficients{t-1} * Mdl.PrivateCovariance * MA.Coefficients{t-1}';
      end

   end
 
else       % Compute a conditional forecast

%
%  Ensure the input model is full-specified.
%
   if any(isnan(Mdl.PrivateConstant))
      error(message('econ:varm:forecast:UnspecifiedConstant'))
   end

   if any(any(isnan([Mdl.AR{:}])))
      error(message('econ:varm:forecast:UnspecifiedAR'))
   end

   if any(isnan(Mdl.PrivateTrend))
      error(message('econ:varm:forecast:UnspecifiedTrend'))
   end

   if any(any(isnan(Mdl.PrivateCovariance)))
      error(message('econ:varm:forecast:UnspecifiedCovariance'))
   end
%
%  Set some flags for later use.
%
   isRegressionIncluded = ~any(strcmpi('X' , parser.UsingDefaults));  % Did the user specify predictors X?
   nPredictors          =  size(X,2);                                 % # of predictor series (X(t) data takes precedence)
   isTrendIncluded      =  any(Mdl.PrivateTrend);
%
%  Now scrub presample responses Y0(t).
%
%  To perform listwise deletion on a multi-variate series, we first convert 
%  the 3-D presample response series Y0(t) to a horizontally concatenated 2-D
%  series. Any row of the combined series with at least one NaN is then 
%  removed. This ensures that the filtered output responses of each path are 
%  the same size and based on the same observation times. 
%
%  However, performing listwise deletion in this manner also implies that in 
%  the presence of missing observations the results obtained from multiple 
%  paths may differ from those obtained from each path individually.
%
   [nRows,~,nPages] = size(Y0);
   Y0               = reshape(Y0, [nRows  (K*nPages)]);

   if any(isnan(Y0(:)))
      Y0 = internal.econ.LagIndexableTimeSeries.listwiseDelete(Y0);
   end
%
%  Determine the origin of the linear time trend to enforce the correct
%  transition from the pre-sample past into the filtered/simulated future.
%
   if isTrendIncluded
      tOrigin = size(Y0,1) - P;
   end

   try
      Y0 = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(P,K*nPages), 'Y0', Y0, P);
   catch exception
      exception.throwAsCaller();
   end

   Y0 = reshape(Y0,[P K nPages]);
%
%  Scrub predictor forecast X(t).
%
   if isRegressionIncluded
 
      if any(isnan(X(:)))    % Check X(t) for missing observations
%
%        When performing a conditional forecast, the rows of predictor forecasts 
%        X(t) and response forecasts YF(t) are associated with the same observation 
%        times, and any rows of YF associated with any NaNs in X are removed
%        via listwise deletion (but NOT vice versa!).
%
         iNonMissing = ~isnan(sum(X,2));
         minimum     =  min([size(X,1) size(Y,1) horizon]);
         Y           =  Y([iNonMissing((1:minimum)) ; true(size(Y,1) - minimum,1)],:,:);
         X           =  internal.econ.LagIndexableTimeSeries.listwiseDelete(X);
      end

      try
%
%        The following code segment ensures we strip data from the end of
%        the X data forecast (i.e., we strip the most distant forecasts).
%
         X = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(horizon,size(X,2)), 'X', X(end:-1:1,:), horizon);
         X = X(end:-1:1,:);
      catch exception
         if strcmp(exception.identifier, 'econ:LagIndexableTimeSeries:checkPresampleData:InvalidPreSampleRows')
            error(message('econ:varm:forecast:InconsistentX', horizon))
         else
            exception.throwAsCaller();
         end
      end
   
      if isempty(Mdl.PrivateBeta) || any(any(isnan(Mdl.PrivateBeta)))
         error(message('econ:varm:forecast:UnspecifiedBeta'))
      else
         if nPredictors ~= size(Mdl.PrivateBeta,2)
            error(message('econ:varm:forecast:InconsistentRegressionSpecification', nPredictors, size(Mdl.PrivateBeta,2)));
         end
      end
   end
%
%  Scrub input response forecasts Y(t).
%
   if size(Y,1) < horizon
      error(message('econ:varm:forecast:InsufficientYFRows', horizon))
   end
%
%  Determine the number of output paths.
%
   if (size(Y0,3) > 1) && (size(Y,3) > 1)
      nPaths = min(size(Y0,3), size(Y,3));
   else
      nPaths = max(size(Y0,3), size(Y,3));
   end
%
%  Now re-size Y0(t) and Y(t).
%
   Y0 = Y0(:,       :,1:min(size(Y0,3),nPaths));
   Y  = Y(1:horizon,:,1:min(size(Y ,3),nPaths));
%
%  Handle the special case of a fully-determined input future response 
%  series (i.e., Y(t) has no NaNs) by simply returning the input series Y(t).
%
   if all(~isnan(Y(:)))
 
      if nargout > 1
         YMSE = repmat({zeros(K)}, horizon, 1);    % Then it's all 0's
      end

   else
%
%     Format the regression and/or trend information into a cell vector as
%     required by SSM.
%      
      if isRegressionIncluded || isTrendIncluded
         fixedTerm = cell(horizon,1);
         for t = 1:horizon
             fixedTerm{t} = Mdl.PrivateConstant;
             if isRegressionIncluded
                fixedTerm{t} = fixedTerm{t} + Mdl.PrivateBeta * X(t,:)';
             end
             if isTrendIncluded
                fixedTerm{t} = fixedTerm{t} + Mdl.PrivateTrend * (t + tOrigin);
             end
         end
      else
         fixedTerm = Mdl.PrivateConstant;
      end
%
%     Now compute conditional forecast and forecast covariances.
%
      if nargout > 1
         [Y,YMSE] = getConditionalForecast(Mdl, Y0, Y, fixedTerm);
      else
         Y = getConditionalForecast(Mdl, Y0, Y, fixedTerm);
      end

   end

end

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function [YF,YMSE] = getConditionalForecast(Mdl, Y0, Y, FixedTerm)
% Conditional forecasts and mean square errors of VAR model responses.
%
% Syntax:
%
%   [Y,YMSE] = getConditionalForecast(Mdl,Y0,Y,FixedTerm)
%
% Input Arguments:
%
%   Mdl - VAR model whose responses Y(t) are forecasted.
%
%   Y0 - Presample response array providing initial values for the forecast. 
%     Y0 is a 2-D matrix (single path) or 3-D (multiple paths) array with 
%     P rows and K columns.
%
%   Y - Time series of forecasted (future) responses to support conditional 
%     forecasting, capturing the state or knowledge of responses as they 
%     evolve from the presample past (Y0) into the future (Y). Although the 
%     number of rows in Y (numperiods) indicates the forecast horizon and 
%     the number of columns must be K = numseries, Y may be a 2-D matrix (single
%     path) or 3-D array (multiple paths) indicating a variable number of 
%     forecast paths.
%
%   FixedTerm - numseries-by-1 constant vector (Mdl.Constant) or numperiods-by-1
%     cell vector of time-varying deterministic terms to capture the sum of
%     the constant, regression component, and linear time trend:
%
%     Mdl.Constant + Mdl.Beta * X(t,:)' + Mdl.Trend * t
%
% Output Arguments:
%
%   Y - numperiods-by-numseries-by-numpaths matrix of minimum mean square 
%     error (MMSE) response forecasts. 
%
%   YMSE - numperiods-by-1 cell vector of mean square errors (MSE) of the
%     response forecasts in Y.

%
% Determine basic array dimensions.
%

P           = Mdl.P;           % P = size(Y0,1) at this point
nY0         = P;               % But store it now because when P = 0 it's re-set to 1 below!
[horizon,K] = size(Y(:,:,1));
nPagesY0    = size(Y0,3);
nPagesY     = size(Y ,3);
nPaths      = max(nPagesY0, nPagesY);

if P == 0
   AR = zeros(K);    % Replicate a VAR(0) model by a zero-valued VAR(1)
   P  = 1;
else
   AR = [Mdl.AR{:}]; % VAR coefficient concatenation
end

%
% Express the K-D VARX(P) model in state-space representation:
%
% For K-by-1 response y(t) and innovation e(t) vectors and an N-by-1 vector
% of predictors x(t), a VARX(P) model has the general form:
%
%   y(t) = c + A1*y(t-1) + A2*y(t-2) + ... + Ap*y(t-P) + B*x(t) + T*t + e(t)
%        = A1*y(t-1) + A2*y(t-2)... + Ap*y(t-P) + F(t) + e(t)
%
% where 
%
%   F(t) = c + B*x(t) + T*t
%
% is a K-by-1 vector of (possibly) time-varying coefficients associated 
% with the remaining non-stochastic terms in the model.
%   
% Now re-write the VARX(P) model as a stacked VARX(1) model and express it
% in state-space representation (see Lutkepohl, pages 612 & 615):
%
% State equation:       Y(t) = A(t) * Y(t-1) + B * u(t)
% Observation equation: y(t) = C * Y(t)
%
% where u(t) is a K-by-1 uncorrelated, unit-variance, Gaussian white noise 
% vector process such that e(t) = L*u(t) and
%
%   Y(t) = [y(t) y(t-1) ... y(t-P+1) 1]'       (K*P + 1)-by-1
%
%          [A1 A2   ...  Ap-1 Ap F(t)           K-by-(K*P + 1)
%   A(t) =   I  0   ...   0    0  0             K-by-(K*P + 1)
%            0  I   ...   0    0  0             K-by-(K*P + 1)
%            .  0   .     0    0  0    =        K-by-(K*P + 1)
%            .  .     .        .  .                 .
%            .  .       .      .  .                 .
%            0  .         I    0  0             K-by-(K*P + 1)
%            0  0   ...   0    0  1 ]           1-by-(K*P + 1)
%                                       ----------------------
%                                       (K*P + 1)-by-(K*P + 1)
%
%   B = [L 0 ... 0 0]'                          (K*P + 1)-by-K
%
%   C = [I 0 ... 0 0]                           K-by-(K*P + 1)
%
% and where L is the K-by-K lower Cholesky factor of the innovations e(t) 
% covariance matrix (Mdl.Covariance).
%

%
% Create the state transition matrix A.
%

A2 = [eye(K*(P-1))  zeros(K*(P-1),K+1)];
A3 = [zeros(1,K*P)  1                 ];

if iscell(FixedTerm)
%
%  Format a time-varying state transition matrix A(t).
%
   A = cell(horizon,1);

   for t = 1:horizon        
       A1   = [AR   FixedTerm{t}(:)];        
       A{t} = [A1 ; A2 ; A3];
   end
else
% 
%  Format a time-invariant state transition matrix A.
%
   A1 = [AR   FixedTerm(:)];
   A  = [A1 ; A2 ; A3];
end

%
% Create the state disturbance loading matrix B.
%

covariance = 0.5 * (Mdl.PrivateCovariance + Mdl.PrivateCovariance');
B1         = cholcov(covariance)';
B2         = zeros(K*(P-1)+1,size(B1,2));
B          = [B1 ; B2];

if size(B1,2) == 0
   warning(message('econ:varm:forecast:DegenerateCovarianceMatrix'))
end

%
% Create the observation sensitivity matrix C.
%

C = [eye(K)  zeros(K,K*(P-1)+1)];

%
% Obtain the initial state distribution from Y0.
%

if nY0 > 0        % Most recent observations prior
   Y0    = permute(Y0(P:-1:P-P+1,:,:), [2 1 3]);
   Cov0  = zeros(K*P+1);
end

%
% Create the state-space model (SSM).
%

if nY0 == 0       % Stationary prior
   SSM = ssm(A, B, C, 'Mean0', [], 'Cov0', []);
else              % Most recent observations prior
   if nPagesY0 == 1
      SSM = ssm(A, B, C, 'Mean0', [Y0(:) ; 1], 'Cov0', Cov0);
   else
      SSM = ssm(A, B, C, 'Cov0', Cov0);
   end
end

%
% Compute conditional forecasts by a Kalman filter (see SSM/FILTER). 
%
% Note that the Kalman smoother (see SSM/SMOOTH) could also be used to
% generate forecasts and covariances under different assumptions.
%

YF = zeros(horizon,K,nPaths);   % Pre-allocate the output forecast array.

for iPath = 1:nPaths
    if (nY0 > 0) && (nPagesY0 > 1)
       YY0       = Y0(:,:,iPath);
       SSM.Mean0 = [YY0(:) ; 1];
    end
    if nPagesY > 1
       [States,~,Output] = filter(SSM, Y(:,:,iPath));
%      [States,~,Output] = smooth(SSM, Y(:,:,iPath));    % Substitute SMOOTH for FILTER
    else
       [States,~,Output] = filter(SSM, Y);
%      [States,~,Output] = smooth(SSM, Y);               % Substitute SMOOTH for FILTER
    end
    YF(:,:,iPath) = States(:,1:K);
end

%
% MSE of the forecast.
%

if nargout > 1
    YMSE = cell(horizon,1);
    for t = 1:horizon
        YMSE{t} = Output(t).FilteredStatesCov(1:K,1:K);
%       YMSE{t} = Output(t).SmoothedStatesCov(1:K,1:K);  % When using the SMOOTH method 
    end
end

end