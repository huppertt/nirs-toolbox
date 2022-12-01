function [E,logL] = infer(Mdl, Y, varargin)
%INFER Infer vector autoregression (VAR) model innovations
%
% Syntax:
%
%   [E,logL] = infer(Mdl,Y)
%   [E,logL] = infer(Mdl,Y,param1,val1,param2,val2)
%
% Description:
%
%   Infer the innovations of a vector autoregression (VAR) model.
%
% Input Arguments:
%
%   Mdl - VAR model created by the VARM constructor or VARM/ESTIMATE method.
%
%   Y - A numobs-by-numseries-by-numpaths 3-D matrix of responses, and the 
%     continuation of the presample response series Y0 (see below). Each page
%     of Y represents a single path of a numseries-dimensional response series,
%     and for all pages observations across any row occur at the same time. 
%     The last observation of any page is the most recent.
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
%        first numpaths pages are used. The default is to strip the first 
%        Mdl.P observations from the beginning of the responses Y, reducing 
%        the effective sample size. The last observation of any page is the 
%        most recent.
%
%  'X'   Predictor data corresponding to a regression component in the 
%        model. X is a 2-D matrix with numpreds columns representing a 
%        numpreds-dimensional time series of predictors, and the last row 
%        contains the most recent observation. The number of observations 
%        (rows) of X required depends upon the specification of presample 
%        responses Y0 (see above). If presample responses Y0 are also 
%        specified, then the number of observations of X must be at least 
%        numobs. Otherwise, the number of observations of X must be at least 
%        numobs - Mdl.P to account for presample stripping. If X has more 
%        observations than necessary, only the most recent observations are 
%        used. The default is an empty series and the estimated VAR model 
%        has no regression component.
%
% Output Arguments:
%
%   E - The inferred model innovations. If presample responses Y0 are 
%     specified, then E is a numobs-by-numseries-by-numpaths matrix the same 
%     size as Y; otherwise the number of observations (rows) in E is 
%     numobs - Mdl.P to account for presample stripping.
%
%   logL - numpaths element vector of loglikelihood objective function 
%     values associated with the model specification Mdl. Each element of 
%     logL is associated with the corresponding path of Y.
%
% Notes:
%
%   o Missing values, indicated by NaNs, are removed from Y and X by listwise 
%     deletion, reducing the effective sample size. That is, the 3-D 
%     responses Y are first converted to an equivalent 2-D series by 
%     horizontal concatenation, and are then merged with the predictors X to 
%     form a composite series in which the last row of each series is assumed 
%     to occur at the same time. Any row of the composite series with at 
%     least one NaN is then removed. This ensures that the inferred output 
%     innovations of each path are the same size and based on the same 
%     observation times, yet also implies that in the presence of missing 
%     observations the results obtained from multiple paths of Y may 
%     differ from those obtained from each path individually.
%
%     Similarly, missing values in the presample data Y0 are also removed by 
%     listwise deletion (i.e., 3-D to 2-D concatenation).
%
%  o  The inclusion of a regression component is based on the presence of 
%     the predictor matrix X. Although each page of the output innovation E 
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
%     Therefore, the times used in the trend component of the output 
%     innovations E(t) are based on times t = t0 + 1, t0 + 2, ..., t0 + numobs.
%     This convention is consistent with the default behavior of model 
%     estimation in which the first Mdl.P responses are stripped, reducing 
%     the effective sample size. If Y0 is unspecified, t0 = 0.
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
% See also ESTIMATE.

% Copyright 2019 The MathWorks, Inc.   

if numel(Mdl) > 1
   error(message('econ:varm:infer:NonScalarModel'))
end

%
% Ensure the input model is full-specified.
%

if any(isnan(Mdl.PrivateConstant))
   error(message('econ:varm:infer:UnspecifiedConstant'))
end

if any(any(isnan([Mdl.AR{:}])))
   error(message('econ:varm:infer:UnspecifiedAR'))
end

if any(isnan(Mdl.PrivateTrend))
   error(message('econ:varm:infer:UnspecifiedTrend'))
end

if any(any(isnan(Mdl.PrivateCovariance)))
   error(message('econ:varm:infer:UnspecifiedCovariance'))
end

%
% Initialize some basic parameters and perform error-checking.
%

K    = size(Mdl.PrivateConstant,1);     % # of response series (numseries) in Y(t)
P    = Mdl.P;                           % # of pre-sample responses need for initialization
Lags = Mdl.PrivateARLagOp.Lags;         % Lags of AR polynomial
Lags = Lags(Lags > 0);                  % Retain only positive lags
AR   = Mdl.AR;                          % Get AR coefficients at positive lags

parser = inputParser;
parser.addRequired ('Y' ,       @(x) validateattributes(x, {'double'}, {'ncols' K 'nonempty'}, '', 'disturbances'));
parser.addParameter('X' ,   [], @(x) validateattributes(x, {'double'}, {'2d'}                , '', 'predictors'));
parser.addParameter('Y0',   [], @(x) validateattributes(x, {'double'}, {'ncols' K           }, '', 'presample responses'));
parser.parse(Y, varargin{:});

Y  = parser.Results.Y;
X  = parser.Results.X;
Y0 = parser.Results.Y0;

%
% Set some flags for later use.
%

isY0Specified        = ~any(strcmpi('Y0', parser.UsingDefaults));   % Did the user specify presample responses Y0?
isRegressionIncluded = ~any(strcmpi('X' , parser.UsingDefaults));   % Did the user specify predictors X?
nPredictors          =  size(X,2);                                  % # of predictor series (X(t) data takes precedence)
isTrendIncluded      =  any(Mdl.Trend);                             % Does the model include a time trend?

%
% Listwise-delete any missing observations from the responses Y(t) and 
% predictors X(t). 
%
% To perform listwise deletion on a multi-variate series, we first convert 
% the 3-D response series Y(t) to a horizontally concatenated 2-D series, 
% then merge Y(t) with the predictors X(t) to form a composite series in which 
% the last row of all paths are assumed to occur at the same time. Any row 
% of the combined series with at least one NaN is then removed. This ensures 
% that the filtered output responses of each path are the same size and 
% based on the same observation times. 
%
% However, performing listwise deletion in this manner also implies that in 
% the presence of missing observations the results obtained from multiple 
% paths may differ from those obtained from each path individually.
%

if any(isnan(Y(:))) || any(isnan(X(:))) 
   [nRows,~,nPages] = size(Y);
   Y                = reshape(Y, [nRows  (K*nPages)]);
   [Y,X]            = internal.econ.LagIndexableTimeSeries.listwiseDelete(Y,X);
   Y                = reshape(Y,[size(Y,1) K nPages]);
end

%
% Now that the input responses have been scrubbed, derive the number of paths 
% and re-format Y(t) such that time corresponds to columns rather than rows
% to allow sequential access to all series at any given time and improved 
% runtime performance.
%

nPaths = size(Y,3);             % # of paths
Y      = permute(Y,[2 1 3]);    % Ensure time corresponds to columns rather than rows

%
% Now scrub any presample responses Y0(t) specified by the user.
%
% When Y0(t) is specified, the following code segment removes missing 
% observations (NaNs) via listwise deletion in the same manner as above for 
% the responses Y(t) (i.e., 3-D to 2-D concatenation), and then ensures 
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
         error(message('econ:varm:infer:InconsistentY0', nPaths))
      else
         exception.throwAsCaller();
      end
   end
   Y0 = reshape(Y0,[P K nPaths]);
   Y0 = permute(Y0,[2 1 3]);           % Ensure time corresponds to columns rather than rows
   T  = size(Y,2);                     % Effective sample size

else

%
%  Since no pre-sample responses are specified, assume the linear trend
%  times are t = 1, 2, ..., T as for model estimation.
%
   if isTrendIncluded
      tOrigin = 0;
   end

   T = size(Y,2) - P;                  % Effective sample size to account for auto-stripping of presample responses
   
end

%
% Missing observations (NaNs) in responses Y(t) and predictors X(t) have 
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
         error(message('econ:varm:infer:InconsistentX', T))
      else
         exception.throwAsCaller();
      end
   end
   X = X';     % Ensure time corresponds to columns rather than rows
   
   if isempty(Mdl.PrivateBeta) || any(any(isnan(Mdl.PrivateBeta)))
      error(message('econ:varm:infer:UnspecifiedBeta'))
   else
      if nPredictors ~= size(Mdl.PrivateBeta,2)
         error(message('econ:varm:infer:InconsistentRegressionSpecification', nPredictors, size(Mdl.PrivateBeta,2)))
      end
   end
   beta = Mdl.PrivateBeta;
end

%
% Extract the linear time trend if necessary.
%

if isTrendIncluded
   trend = Mdl.PrivateTrend;
end

%
% Prepend the presample responses and pre-allocate the innovations such 
% that E(t) = Y(t) - C.
%

if isY0Specified
   Y = [Y0  Y];       % Prepend presample responses for convenience
end

E = Y - Mdl.PrivateConstant;

%
% Now infer the innovations E(t) from the responses Y(t). 
%
% Note that at this point the innovations E(t) already account for the current 
% response Y(t) and constant C.
%

for iPath = 1:nPaths

    for t = (P + 1):(P + T)
% 
%       Subtract autoregressive terms from the innovations.
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

	   end
end

E = permute(E(:,P+1:end,:), [2 1 3]);  % Convert back to time series format

if nargout > 1
   logL = getLogL(Mdl, E);             % Compute the Gaussian log-likelihood 
end

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function logL = getLogL(Mdl, E)
%GETLOGL Compute the log-likelihood of a multivariate Gaussian objective
%
% Syntax:
%
%   logL = getLogL(Mdl,E)
%
% Inputs:
%
%   Mdl - VAR model, as produced by VARM constructor or VARM/ESTIMATE.
%
%   E - The inferred model innovations. E is a numobs-by-numseries-by-numpaths
%     matrix of residuals.
%
% Output:
%
%   logL - numpaths element vector of loglikelihood objective function 
%     values associated with the model specification Mdl. Each element of 
%     logL is associated with the corresponding path (i.e., page) of E.
%
% Note:
%
% This function is taken directly from the relevant segments of the undocumented 
% function STATMVNROBJ in the Statistics Toolbox. For more details, see 
% \matlab\toolbox\stats\stats\private\statmvnrobj.m.
%

[T,K,nPaths] = size(E);
R            = chol(Mdl.Covariance);     % Upper triangular Cholesky factor 
LogTwoPi     = log(2.0 * pi);
LogDetCovar  = 2.0 * sum(log(diag(R)));  % log(det(Mdl.Covariance))
logL         = nan(1,nPaths);

for iPath = 1:nPaths
    logL(iPath) = -0.5 * sum(sum((E(:,:,iPath) / R).^2,2));
    logL(iPath) =  logL(iPath) - 0.5 * T * (K * LogTwoPi + LogDetCovar);
end

end