function [Y,E] = simulate(Mdl, numobs, varargin)
%SIMULATE Simulate vector autoregression (VAR) model responses
%
% Syntax:
%
%   [Y,E] = simulate(Mdl,numobs)
%   [Y,E] = simulate(Mdl,numobs,param1,val1,param2,val2,...)
%
% Description:
%
%   Simulate sample paths of responses and innovations of a vector 
%   autoregression (VAR) model.
%
% Input Arguments:
%
%   Mdl - VAR model created by the VARM constructor or VARM/ESTIMATE method.
%
%   numobs - Positive integer indicating the number of observations (rows)
%     simulated for each path of the output responses (Y) and innovations (E).
%
% Optional Input Parameter Name/Value Pairs:
%
%'NumPaths'  Positive integer indicating the number of sample paths (pages) 
%        generated for all simulated outputs. The default is 1.
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
%        X must equal or exceed the number of observations simulated (numobs).
%        When the number of observations in X exceeds the number necessary, 
%        only the most recent observations are used. By default the simulation
%        includes no regression component.
%
%  'YF'  Time series of future responses to support conditional simulation. 
%        YF is a 2-D or 3-D matrix with numseries columns, and captures the 
%        state or knowledge of responses as they evolve from the presample 
%        past (Y0) into the future. The first row contains the responses 
%        one period into the future, the second row two periods into the 
%        future, and so forth. As a 2-D matrix, known values in YF are treated
%        as deterministic and are applied to each of numpaths output paths. 
%        As a 3-D matrix, YF provides deterministic values for each output 
%        path and must have at least numpaths pages. YF may have any number 
%        of rows, provided at least numobs rows exist to cover the simulation
%        horizon. If the number of rows exceeds numobs, then only the first 
%        numobs rows are used. If the number of pages exceeds numpaths, then 
%        only the first numpaths pages are used. By default, YF is an array 
%        of NaN values indicating a complete lack of knowledge of the future 
%        state of all simulated responses, and the output responses Y (see 
%        below) are obtained from a conventional Monte Carlo simulation.
%
% Output Arguments:
%
%   Y - numobs-by-numseries-by-numpaths matrix of simulated responses, and 
%     the continuation of the presample series Y0 (see above). If YF is also
%     an input, indicating a conditional simulation, then the output Y will 
%     incorporate the additional information known about future responses.
%
%   E - numobs-by-numseries-by-numpaths matrix of model innovations. 
%
% Notes:
%
%   o Missing data values, indicated by NaNs, are removed from X by listwise 
%     deletion (i.e., any row in X with at least one NaN is removed), reducing
%     the effective sample size. When performing a conditional simulation (see
%     below) there can be no missing values in any of the most recent numobs 
%     observations of X. 
%
%     Similarly, missing values in the presample data Y0 are also removed by 
%     listwise deletion. However, Y0 is generally a 3-D time series array, 
%     and to ensure that the simulated output responses of each path are the 
%     same size and based on the same observation times, Y0 is first converted
%     to an equivalent 2-D series by horizontal concatenation, and any row 
%     of the concatenated series with at least one NaN is removed. This 
%     convention also implies that in the presence of missing observations 
%     the results obtained by simulating multiple paths of may differ from 
%     those obtained by simulating each path individually.
%
%   o When YF is specified as an input, NaNs also indicate missing values, 
%     but also indicate the pattern of partial knowledge of the state of 
%     future responses. Non-missing values in YF serve as deterministic future
%     responses known in advance (e.g., set by policy), and in such cases 
%     any corresponding missing NaN values are simulated conditional on
%     non-missing information. 
%
%     When missing information is included in the input responses YF, the 
%     simulation proceeds in three steps. First, at any time t the innovations
%     E(t) are inferred (inverse filtered) from the input future responses 
%     YF(t), and the pattern of missing information (NaNs) that appears in 
%     YF(t) also appears in E(t). Given the missingness pattern found in E(t), 
%     unknown elements of E(t) are then imputed by sampling from a Gaussian 
%     distribution conditional on the known elements of E(t). Given the 
%     updated innovations, the unknown responses are then updated. These 
%     three steps are repeated until a complete output response series Y(t) 
%     is obtained.
%
%     When YF is unspecified, or otherwise contains completely missing values,
%     the output innovations E(t) and responses Y(t) are obtained from a
%     conventional (unconditional) Monte Carlo simulation. Given an 
%     independent, numseries-dimensional standard Gaussian vector Z(t) and 
%     the lower triangular Cholesky factor (L) of the model covariance matrix 
%     (Mdl.Covariance) the innovations covariance structure at any time t 
%     is such that E(t) = L*Z(t).
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
% See also FILTER.

% Copyright 2017 The MathWorks, Inc.

if numel(Mdl) > 1
   error(message('econ:varm:simulate:NonScalarModel'))
end

%
% Initialize some basic parameters and perform error-checking.
%

K = size(Mdl.PrivateConstant,1);   % # of responses (numseries) in Y(t)

parser = inputParser;
parser.addRequired ('RequiredNumObs',      @(x) validateattributes(x, {'double'}, {'scalar' 'integer' '>' 0}, '', 'number of observations'));
parser.addParameter('NumPaths'      ,   1, @(x) validateattributes(x, {'double'}, {'scalar' 'integer' '>' 0}, '', 'number of paths'));
parser.addParameter('X'             ,  [], @(x) validateattributes(x, {'double'}, {'2d'}          , '', 'predictors'));
parser.addParameter('Y0'            ,  [], @(x) validateattributes(x, {'double'}, {'ncols'   K   }, '', 'presample responses'));
parser.addParameter('YF'            , nan, @(x) validateattributes(x, {'double'}, {'ncols'   K   }, '', 'future responses'));
parser.parse(numobs, varargin{:});

T      = parser.Results.RequiredNumObs;    % # of observations = sample size
nPaths = parser.Results.NumPaths;
X      = parser.Results.X;
Y0     = parser.Results.Y0;
Y      = parser.Results.YF;

%
% Determine whether an unconditional or conditional simulation is selected.
%
% A conditional simulation is indicated the presence of an input response 
% series Y(t) whose elements are NOT all NaNs.
%

if (any(strcmpi('YF', parser.UsingDefaults)) || all(isnan(Y(:))))
%
%  The user requested a conventional (unconditional) Monte Carlo simulation, 
%  so just remove optional parameters specific to SIMULATE and call FILTER.
%
   varargin(1:2:end) =  cellstr(varargin(1:2:end));
   indices           = ~(strncmpi('NumPaths', varargin(1:2:end), 3) | strcmpi('YF', varargin(1:2:end)));
   varargin          =  reshape(varargin, 2, numel(varargin)/2);
   varargin          =  varargin(:,indices);
   
   try
%
%     The following simulates additional observations of Z(t) to compensate 
%     for subsequent missing data (if any) removal from X(t) in FILTER. However, 
%     depending upon where the missing observations appear in X(t), this may 
%     not be necessary and so any excess observations are trimmed before 
%     returning.
%
      if nargout > 1
         [Y,E] = filter(Mdl, randn(T + sum(isnan(sum(X,2))),K,nPaths), varargin{:});
         E     = E(1:T,:,:);     % Retain the first T observations if necessary
      else
         Y = filter(Mdl, randn(T + sum(isnan(sum(X,2))),K,nPaths), varargin{:});
      end
   catch exception
      exception.throwAsCaller();
   end
   Y = Y(1:T,:,:);               % Retain the first T observations if necessary
   return
end

%
% Perform a conditional Monte Carlo simulation.
%
% Note that the following code path will work for ANY arbitrary missingness
% pattern, including situations in which the input future response series Y(t) 
% is a mix of NaN (unknown) and non-NaN (known) values, all NaNs (Y(t) is 
% completely unknown), and no NaNs (Y(t) is fully-determined).

%
% Ensure the input model is full-specified.
%

if any(isnan(Mdl.PrivateConstant))
   error(message('econ:varm:simulate:UnspecifiedConstant'))
end

if any(any(isnan([Mdl.AR{:}])))
   error(message('econ:varm:simulate:UnspecifiedAR'))
end

if any(isnan(Mdl.PrivateTrend))
   error(message('econ:varm:simulate:UnspecifiedTrend'))
end

if any(any(isnan(Mdl.PrivateCovariance)))
   error(message('econ:varm:simulate:UnspecifiedCovariance'))
end

P    = Mdl.P;                           % # of pre-sample responses needed for initialization
Lags = Mdl.PrivateARLagOp.Lags;         % Lags of AR polynomial
Lags = Lags(Lags > 0);                  % Retain only positive lags
AR   = Mdl.AR;                          % Get AR coefficients at positive lags

%
% Set some flags for later use.
%

isY0Specified        = ~any(strcmpi('Y0', parser.UsingDefaults));   % Did the user specify presample responses Y0?
isRegressionIncluded = ~any(strcmpi('X' , parser.UsingDefaults));   % Did the user specify predictors X?
nPredictors          =  size(X,2);                                  % # of predictor series (X(t) data takes precedence)
isTrendIncluded      =  any(Mdl.PrivateTrend);       

%
% Now scrub any presample responses Y0(t) specified by the user, or
% auto-generate accordingly.
%
% When Y0(t) is specified, the following code segment removes missing 
% observations (NaNs) via listwise deletion by converting the input 3-D series 
% to a horizontally concatenated 2-D series. Any row of the combined series 
% with at least one NaN is then removed. This ensures that the output responses 
% of each path are the same size and based on the same observation times. 
%
% Y0(t) is also checked to ensure it is of sufficient length (P).
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
%  and listwise-delete any rows with at least one NaN. 
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
         error(message('econ:varm:simulate:InconsistentY0',nPaths))
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
% Scrub any predictors X(t) indicating the inclusion of a regression component.
%
% The following code segment removes missing observations (NaNs) in predictors 
% X(t) by listwise deletion, then ensures that the predictors X(t) are of 
% sufficient length. The following code segment also ensures consistency
% between the predictor coefficients of the model and the number of predictors 
% in X(t).
%
% After X(t) is scrubbed, the data is re-formatted such that time corresponds 
% to columns rather than rows. 
%

if isRegressionIncluded
 
   if any(isnan(X(:)))    % Check X(t) for missing observations
%
%     When performing a conditional simulation with a regression component, 
%     recall that the FIRST T observations of future responses YF(t) are
%     used, whereas the LAST T (most recent) observations of predictors X(t) 
%     are used.
%
%     To maintain consistent synchronization of observations times between 
%     future responses YF(t) and predictors X(t), any missing observations in 
%     X MUST occur at the beginning of X(t) and there cannot be any NaNs in 
%     the most recent (the last) T rows of X.
%
      if sum(~isnan(sum(X((end - min([size(X,1) size(Y,1) T]) + 1):end,:),2))) < T
         error(message('econ:varm:simulate:InvalidMissingnessPatternX', T))
      end
      X = internal.econ.LagIndexableTimeSeries.listwiseDelete(X);
   end

   try 
      X = internal.econ.LagIndexableTimeSeries.checkPresampleData(zeros(T,size(X,2)), 'X', X, T);
   catch exception
      if strcmp(exception.identifier, 'econ:LagIndexableTimeSeries:checkPresampleData:InvalidPreSampleRows')
         error(message('econ:varm:simulate:InconsistentX',T))
      else
         exception.throwAsCaller();
      end
   end    
   X = X';     % Ensure time corresponds to columns rather than rows
   
   if isempty(Mdl.PrivateBeta) || any(any(isnan(Mdl.PrivateBeta)))
      error(message('econ:varm:simulate:UnspecifiedBeta'))
   else
      if nPredictors ~= size(Mdl.PrivateBeta,2)
         error(message('econ:varm:simulate:InconsistentRegressionSpecification', nPredictors, size(Mdl.PrivateBeta,2)))
      end
   end
   beta = Mdl.PrivateBeta;
end

%
% Scrub the input future responses Y(t).
%

if size(Y,1) < T
   error(message('econ:varm:simulate:InsufficientYFRows', T))
end
   
Y = Y(1:T,:,:);                   % Re-size the # of observations

if (size(Y,3) == 1) && (nPaths > 1)
   Y = repmat(Y, [1 1 nPaths]);   % Re-size the # of paths
end
if size(Y,3) < nPaths
   error(message('econ:varm:simulate:InsufficientYFPaths', nPaths))
else
   Y = Y(:,:,1:nPaths);           % Re-size the # of paths
end

Y = permute(Y,[2 1 3]);           % Ensure time corresponds to columns rather than rows

%
% Prepend the presample responses Y0(t) to Y(t) and pre-allocate the innovations E(t).
%

Y = [Y0  Y];
E = zeros(size(Y));

%
% Extract the linear time trend if necessary.
%

if isTrendIncluded
   trend = Mdl.PrivateTrend;
end

%
% Now infer (inverse-filter) the innovations E(t) from the responses Y(t), 
% impute any missing values found in E(t), and then simulate (filter) the
% innovations E(t) to update the responses Y(t).
%

constant = Mdl.PrivateConstant;
sigma    = Mdl.PrivateCovariance;
I        = eye(K);

for iPath = 1:nPaths

    for t = (P + 1):(P + T)
% 
%       Subtract the constant from the current response.
%  
        E(:,t,iPath) = Y(:,t,iPath) - constant;
% 
%       Subtract autoregressive terms.
%
     			for Lag = Lags
					       E(:,t,iPath) = E(:,t,iPath) - AR{Lag} * Y(:,t - Lag,iPath);
        end
%
%       Subtract the regression component.
%
        if isRegressionIncluded
           E(:,t,iPath) = E(:,t,iPath) - beta * X(:,t-P);
        end
%
%       Subtract the linear time trend.
%
        if isTrendIncluded
           E(:,t,iPath) = E(:,t,iPath) - trend * (t - P + tOrigin);
        end
%
%       Generate conditional Gaussian draws if necessary using the Gaussian 
%       Conditioning Formula (see Glasserman page 65, equations 2.24 and 2.25, 
%       or Lutkepohl page 678).
%
%       Note that since the pre-sample responses Y0(t) and predictors X(t) 
%       have no missing (NaN) values, the ONLY way that the current innovation
%       E(t) can have NaN values is if NaNs appear in the current response Y(t).
%
        i1 = isnan(E(:,t,iPath));       % Missing value (NaN) indicator

        if any(i1)
           i2 = ~i1;                    % Non-missing value (non-NaN) indicator
%
%          Compute the conditional covariance of the unknown innovations
%          and sample the unknown E(t)'s from the conditional distribution.
%
           term          = sigma(i1,i2) * (sigma(i2,i2) \ I(i2,i2));
           C             = sigma(i1,i1) - term * sigma(i2,i1);  % Conditional covariance
           L             = chol(C,'lower');                     % Lower Cholesky factor
           E(i1,t,iPath) = L * randn(sum(i1),1) + term * E(i2,t,iPath);
%
%          Update the current response Y(t) to reflect the current innovation E(t).
%
           Y(:,t,iPath) = E(:,t,iPath) + constant;
           for Lag = Lags
					          Y(:,t,iPath) = Y(:,t,iPath) + AR{Lag} * Y(:,t - Lag,iPath);
           end
           if isRegressionIncluded
              Y(:,t,iPath) = Y(:,t,iPath) + beta * X(:,t-P);
           end
           if isTrendIncluded
              Y(:,t,iPath) = Y(:,t,iPath) + trend * (t - P + tOrigin);
           end
        end
        
    end

end

%
% Convert Y(t) and E(t) back to conventional time series format.
%
% For reference, in a large sample experiment in which a subset of future 
% responses are completely unknown while others are completely known, the 
% conditional covariance (C, computed above) should be approximately 
%
% C ~ cov(E(:,i1) - E(:,i2)*term')
%
% while the mean should be mean(E(:,i1) - E(:,i2)*term').
%

Y = permute(Y(:,(P+1):(P+T),:), [2 1 3]);
E = permute(E(:,(P+1):(P+T),:), [2 1 3]);

end