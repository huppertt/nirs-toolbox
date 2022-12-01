function [Y,E] = filter(Mdl,Y,varargin)
%FILTER Filter vector autoregression (VAR) model disturbances
%
% Syntax:
%
%   [Y,E] = filter(Mdl,Z)
%   [Y,E] = filter(Mdl,Z,param1,val1,param2,val2,...)
%
% Description:
%
%   Filter user-specified disturbances to produce responses and innovations
%   of a vector autoregression (VAR) model.
%
% Input Arguments:
%
%   Mdl - VAR model created by the VARM constructor or VARM/ESTIMATE method.
%
%   Z - A numobs-by-numseries-by-numpaths 3-D matrix of disturbances associated 
%     with the model innovations process E(t). By default, the optional input
%     flag 'Scale' (see below) is true, and the lower triangular Cholesky 
%     factor (L) of the model covariance matrix (Mdl.Covariance) is used to 
%     scale the innovations such that E(t) = L*Z(t). When the 'Scale' flag is 
%     false, no scaling occurs and E(t) = Z(t). Each page of Z represents a 
%     single path of the underlying disturbance series, and for all pages 
%     observations across any row occur at the same time. The last observation
%     of any page is the most recent.
%
% Optional Input Parameter Name/Value Pairs:
%
%  'Y0'  Presample response data providing initial values for the model. Y0 
%        is a 2-D or 3-D matrix with numseries columns. As a 2-D matrix, 
%        initial values in Y0 are applied to each output path and all paths 
%        evolve from common initial conditions. As a 3-D matrix, Y0 provides 
%        initial conditions for each output path and must have at least 
%        numpaths pages. Y0 may have any number of rows, provided at least 
%        Mdl.P observations exist to initialize the model. If the number of 
%        rows exceeds Mdl.P, then only the most recent Mdl.P observations 
%        are used. If the number of pages exceeds numpaths, then only the 
%        first numpaths pages are used. By default, any necessary presample 
%        observations are set to the unconditional mean for stationary VAR 
%        processes, and to zero if the process is non-stationary or contains 
%        a regression component. The last observation of any page is the most 
%        recent.
%
%  'X'   Predictor data corresponding to a regression component in the 
%        model. X is a 2-D matrix with numpreds columns representing a 
%        numpreds-dimensional time series of predictors, and the last row 
%        contains the most recent observation. The number of observations in 
%        X must equal or exceed the number of observations in Z (numobs). 
%        When the number of observations in X exceeds the number necessary, 
%        only the most recent observations are used. The default is an empty 
%        series and the model has no regression component.
%
%'Scale' Scalar logical flag to indicate whether the input disturbances Z(t)
%        are scaled by the lower triangular Cholesky factor (L) of the model 
%        covariance matrix (Mdl.Covariance). If true, the input disturbances
%        Z(t) are scaled to produce the output innovations, E(t) = L*Z(t).
%        If false, the input disturbances are simply assigned to the output 
%        innovations, E(t) = Z(t). The default is true.
%
% Output Arguments:
%
%   Y - numobs-by-numseries-by-numpaths matrix of filtered responses, and the
%     continuation of the presample series Y0 (see above).
%
%   E - numobs-by-numseries-by-numpaths matrix of model innovations. If 'Scale'
%     is true then E(t) = L*Z(t), otherwise E(t) = Z(t).
%
% Notes:
%
%   o The FILTER function generalizes the SIMULATE function. Both functions 
%     filter a series of disturbances to produce output responses Y(t) and 
%     innovations E(t). However, whereas SIMULATE generates a series of
%     mean-zero, unit-variance, independent Gaussian disturbances Z(t) to 
%     form innovations E(t) = L*Z(t), FILTER allows users to directly specify 
%     their own disturbances.
%
%   o Missing values, indicated by NaNs, are removed from Z and X by listwise 
%     deletion, reducing the effective sample size. That is, the 3-D 
%     disturbances Z are first converted to an equivalent 2-D series by 
%     horizontal concatenation, and are then merged with the predictors X to 
%     form a composite series in which the last row of each series is assumed 
%     to occur at the same time. Any row of the composite series containing
%     at least one NaN is then removed. This ensures that the filtered output
%     responses of each path are the same size and based on the same 
%     observation times, yet also implies that in the presence of missing 
%     observations, the results obtained by filtering multiple paths of Z 
%     may differ from those obtained by filtering each path individually.
%
%     Similarly, missing values in the presample data Y0 are also removed by 
%     listwise deletion (i.e., 3-D to 2-D concatenation).
%
%  o  The inclusion of a regression component is based on the presence of 
%     the predictor matrix X. Although each page of the output response Y 
%     represents a different path, the predictor matrix X represents a 
%     single path of a multivariate time series. When the model includes a 
%     regression component, the entire predictor matrix X is applied to 
%     every path of the output response series Y.
%
%  o  Although only the first Mdl.P observations of presample responses Y0(t)
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
%     sample size. If Y0 is unspecified, t0 = 0.
%
% Example:
%
%   o Simulate-to-Filter Round-Trip:
%
%     Consider a "round-trip" example to illustrate the relationship between
%     the simulate and filter functions.
%
%     First, fit a VAR(2) model to the 4-D Danish data obtained from [2]. 
%     Then, simulate a single path of the fitted model using the historical 
%     data as initial values. Compare the output of simulate to that of
%     filter:
%
%     load Data_JDanish
%     Mdl = estimate(varm(4,2), Data);
%
%     rng default
%     Y1 = simulate(Mdl, 100, 'Y0', Data);
%
%     rng default
%     Z = randn(100,4);
%     Y2 = filter(Mdl, Z, 'Y0', Data);
%
%     In the above round-trip example, the outputs of SIMULATE and FILTER 
%     are identical.
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

% Copyright 2018 The MathWorks, Inc.

if numel(Mdl) > 1
   error(message('econ:varm:filter:NonScalarModel'))
end

%
% Ensure the input model is full-specified.
%

if any(isnan(Mdl.PrivateConstant))
   error(message('econ:varm:filter:UnspecifiedConstant'))
end

AR = Mdl.AR;                            % Get AR coefficients at positive lags

if any(any(isnan([AR{:}])))
   error(message('econ:varm:filter:UnspecifiedAR'))
end

if any(isnan(Mdl.PrivateTrend))
   error(message('econ:varm:filter:UnspecifiedTrend'))
end

if any(any(isnan(Mdl.PrivateCovariance)))
   error(message('econ:varm:filter:UnspecifiedCovariance'))
end

%
% Initialize some basic parameters and perform error-checking.
%
% Note on Memory Management:
%
% Although the input disturbances are documented as "Z", they are immediately 
% assigned to the output response variable "Y" in an effort to avoid creating
% unnecessary arrays. Therefore, in what follows the variable "Z" does not 
% appear below, only "Y".
%

K    = size(Mdl.PrivateConstant,1);     % # of response series (numseries) in Y(t)
P    = Mdl.P;                           % # of pre-sample responses need for initialization
Lags = Mdl.PrivateARLagOp.Lags;         % Lags of AR polynomial
Lags = Lags(Lags > 0);                  % Retain only positive lags

parser = inputParser;
parser.addRequired ('Z'    ,       @(x) validateattributes(x, {'double'}          , {'ncols' K 'nonempty'}, '', 'disturbances'));
parser.addParameter('X'    ,   [], @(x) validateattributes(x, {'double'}          , {'2d'}                , '', 'predictors'));
parser.addParameter('Y0'   ,   [], @(x) validateattributes(x, {'double'}          , {'ncols'       K     }, '', 'presample responses'));
parser.addParameter('Scale', true, @(x) validateattributes(x, {'double' 'logical'}, {'scalar'   'integer'}, '', 'scale flag'));

parser.parse(Y, varargin{:});

Y        = parser.Results.Z;            % Assign Y(t) = Z(t)
X        = parser.Results.X;
Y0       = parser.Results.Y0;
isScaled = logical(parser.Results.Scale);

%
% Listwise-delete any missing observations from the disturbances Z(t) and 
% predictors X(t). 
%
% To perform listwise deletion on a multi-variate series, we first convert 
% the 3-D disturbance series Z(t) to a horizontally concatenated 2-D series, 
% then merge Z(t) with the predictors X(t) to form a composite series in which 
% the last row of all paths are assumed to occur at the same time. Any row 
% of the combined series with at least one NaN is then removed. This ensures 
% that the filtered output responses of each path are the same size and 
% based on the same observation times. 
%
% However, performing listwise deletion in this manner also implies that in 
% the presence of missing observations the results obtained by filtering
% multiple paths may differ from those obtained by filtering each path
% individually.
%

if any(isnan(Y(:))) || any(isnan(X(:)))    % Checking Z(t) and X(t)!
   [nRows,~,nPages] = size(Y);
   Y                = reshape(Y, [nRows  (K*nPages)]);
   [Y,X]            = internal.econ.LagIndexableTimeSeries.listwiseDelete(Y,X);
   Y                = reshape(Y,[size(Y,1) K nPages]);
end

%
% Now that the input disturbances have been scrubbed, derive the effective 
% sample size and the number of filtered paths, and re-format such that time 
% corresponds to columns rather than rows to allow sequential access to 
% all series at any given time and improved runtime performance.
%

T      = size(Y,1);             % Effective sample size
nPaths = size(Y,3);             % # of paths
Y      = permute(Y,[2 1 3]);    % Ensure time corresponds to columns rather than rows
   
%
% Set some flags for later use.
%

isY0Specified        = ~any(strcmpi('Y0', parser.UsingDefaults));   % Did the user specify presample responses Y0?
isRegressionIncluded = ~any(strcmpi('X' , parser.UsingDefaults));   % Did the user specify predictors X?
nPredictors          =  size(X,2);                                  % # of predictor series (X(t) data takes precedence)
isTrendIncluded      =  any(Mdl.Trend);                             % Does the model include a time trend?

%
% Missing observations (NaNs) in disturbances Z(t) and predictors X(t) have 
% been removed by listwise deletion, so now ensure that the predictors X(t)
% are of sufficient length. The following code segment also ensures consistency
% between the predictor coefficients of the model and the number of predictors 
% in X(t).
%
% After the time series data is scrubbed, the data is re-formatted such
% that time corresponds to columns rather than rows. 
%

if isRegressionIncluded
   try 
      X = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(T,size(X,2)), 'X', X, T);
   catch exception
      if strcmp(exception.identifier, 'econ:LagIndexableTimeSeries:checkPresampleData:InvalidPreSampleRows')
         error(message('econ:varm:filter:InconsistentX',T))
      else
         exception.throwAsCaller();
      end
   end    
   X = X';     % Ensure time corresponds to columns rather than rows
   
   if isempty(Mdl.PrivateBeta) || any(any(isnan(Mdl.PrivateBeta)))
      error(message('econ:varm:filter:UnspecifiedBeta'))
   else
      if nPredictors ~= size(Mdl.PrivateBeta,2)
         error(message('econ:varm:filter:InconsistentRegressionSpecification', nPredictors, size(Mdl.PrivateBeta,2)))
      end
   end
   beta = Mdl.PrivateBeta;
end

%
% Now scrub any presample responses Y0(t) specified by the user, or
% auto-generate accordingly.
%
% When Y0(t) is specified, the following code segment removes missing 
% observations (NaNs) via listwise deletion in the same manner as above for 
% the disturbances Z(t) (i.e., 3-D to 2-D concatenation), and then ensures 
% the presample responses Y0(t) are of sufficient length.
%
% After the time series data is scrubbed, the data is re-formatted such that
% time corresponds to columns rather than rows. 
%

if isY0Specified

%
%  Presample responses are specified, so remove any missing observations via
%  listwise deletion. 
%
%  Once again, to perform listwise deletion on multiple paths of a multi-variate 
%  series, convert the 3-D series to a horizontally concatenated 2-D series 
%  and listwise delete any rows with at least one NaN. 
%
   [nRows,~,nPages] = size(Y0);
   Y0               = reshape(Y0, [nRows  (K*nPages)]);
   if any(isnan(Y0(:)))
      Y0 = internal.econ.LagIndexableTimeSeries.listwiseDelete(Y0);
   end
   if nPages == 1
      Y0 = repmat(Y0, 1, nPaths);      % Allow a single presample path to initialize all paths
   end
%
%  Determine the origin of the linear time trend to enforce the correct
%  transition from the pre-sample past into the filtered/simulated future.
%
   if isTrendIncluded
      tOrigin = size(Y0,1) - P;
   end
   try
      Y0 = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(P,K*nPaths), 'Y0', Y0, P);
   catch exception
      if strcmp(exception.identifier, 'econ:LagIndexableTimeSeries:checkPresampleData:InvalidPreSampleMatrix')
         error(message('econ:varm:filter:InconsistentY0',nPaths))
      else
         exception.throwAsCaller();
      end
   end      
   Y0 = reshape(Y0,[P K nPaths]);
   Y0 = permute(Y0,[2 1 3]);           % Ensure time corresponds to columns rather than rows

else

%
%  No explicit presample responses are specified, so assign defaults.
%
%  VAR models with stationary AR polynomials, without linear time trends or 
%  regression components, are assigned the unconditional mean of the stationary
%  process. Otherwise, zeros are assigned.
%
   if isStable(Mdl.PrivateARLagOp) && ~isTrendIncluded && (nPredictors == 0)     
      summation = eye(K);
      for Lag = Lags
          summation = summation - AR{Lag};   % (I - A1 - A2 - ... - Ap)
      end
      try
         Y0 = repmat(summation \ Mdl.PrivateConstant, 1, P, nPaths);
      catch
         Y0 = zeros(K,P,nPaths);
      end
   else
      Y0 = zeros(K,P,nPaths);
   end
%
%  Since no pre-sample responses are specified, assume the linear trend
%  times are t = 1, 2, ..., T as for model estimation.
%
   if isTrendIncluded
      tOrigin = 0;
   end

end

%
% Extract the linear time trend if necessary.
%

if isTrendIncluded
   trend = Mdl.PrivateTrend;
end

%
% Compute the output innovations E(t), but store them in the output responses 
% Y(t) to save memory.
%

if isScaled
%
%  Induce covariance structure on the i.i.d. input disturbances Z(t).
%
   L = chol(Mdl.PrivateCovariance,'lower');   % Lower Cholesky factor

   for iPath = 1:nPaths
       Y(:,:,iPath) = L * Y(:,:,iPath);       % E(t) = R'*Z(t) = L*Z(t)
   end

end

%
% Save the innovations E(t) in time series format such that time corresponds 
% to rows rather than columns (i.e., in conventional MATLAB time series format).
%

if nargout > 1
   E = permute(Y,[2 1 3]);
end

%
% Prepend the presample responses and add the constant to the innovations.
%

Y = [Y0  (Y + Mdl.PrivateConstant)];

%
% Now filter the innovations E(t) to produce the responses Y(t). 
%
% Note that at this point the responses Y(t) already include the current 
% innovation E(t) and constant C.
%

for iPath = 1:nPaths

    for t = (P + 1):(P + T)
% 
%       Add autoregressive terms to the response.
%
     			for Lag = Lags
					       Y(:,t,iPath) = Y(:,t,iPath) + AR{Lag} * Y(:,t - Lag,iPath);
        end
%
%       Include a regression component.
%
        if isRegressionIncluded
           Y(:,t,iPath) = Y(:,t,iPath) + beta * X(:,t-P);
        end
%
%       Include a linear time trend.
%
        if isTrendIncluded
           Y(:,t,iPath) = Y(:,t,iPath) + trend * (t - P + tOrigin);
        end

	   end
end

Y = permute(Y(:,P+1:end,:), [2 1 3]);  % Convert back to time series format

end