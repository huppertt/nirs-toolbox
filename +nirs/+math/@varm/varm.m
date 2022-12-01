classdef (Sealed) varm < matlab.mixin.CustomDisplay
%VARM Create vector autoregression (VAR) models
%
% Syntax:
%
%   Mdl = varm(numseries,numlags)
%   Mdl = varm(name1,value1,name2,value2,...)
%
% Description:
%
%   Create a VAR(p) model by specifying the number of time series and the 
%   number of lagged responses (short-hand syntax), or a list of parameter 
%   name-value pairs (long-hand syntax). 
%
%   Either syntax creates a VAR(p) model of the form
%
%   y(t) = c + B1*y(t-1) + ... + Bp*y(t-p) + D*x(t) + T*t + e(t)
%
%   for responses y(t), predictors x(t), and innovations e(t).
%
% Input Arguments (Short-Hand Syntax):
%
%   numseries - Positive integer number of time series, and the dimensionality 
%       of the responses y(t) and innovations e(t).
%
%   numlags - Nonnegative integer number of lagged responses included in the 
%       VAR(p) model. In any given VAR(p) model, numlags = p.
%
% Input Arguments (Parameter Name/Value Pairs):
%
% 'Constant'    numseries-by-1 vector of constants (intercepts), denoted as 
%               c in the VAR equation. The default is a column vector of 
%               NaNs.
%
% 'AR'          Autoregressive coefficients associated with lagged responses, 
%               denoted as B1, B2, ..., Bp in the VAR equation. When specified
%               without corresponding Lags (see below), AR is a p-element 
%               cell vector of numseries-by-numseries coefficient matrices B1, 
%               B2, ..., Bp at lags 1, 2, ..., p. When specified with Lags 
%               AR is a commensurate length cell vector of coefficients 
%               associated with the lags in Lags. The default is a cell vector
%               of NaNs. 
%
% 'Lags'        Vector of unique, positive, integer lags associated with the 
%               autoregressive cell vector of coefficients. The default is a 
%               vector of integers 1, 2, ..., p the same length as AR. 
%
% 'Trend'       numseries-by-1 vector of linear time trends, denoted as T in
%               the VAR equation. The default is a column vector of zeros
%               (no time trend).
%
% 'Beta'        numseries-by-numpreds regression coefficient matrix associated
%               with numpreds predictors in x(t), denoted as D in the VAR 
%               equation. The default is a numseries-by-0 empty matrix.
%
% 'Covariance'  numseries-by-numseries positive definite covariance matrix 
%               of the innovations e(t). The default is a matrix of NaNs.
%
% 'Description' Description of the model, specified as a string. The default 
%               is a summary of the parametric form of the model.
%
% 'SeriesNames' String vector of response series names with numseries elements.
%               The default is a string vector "Y1", "Y2", "Y3", and so forth.
%
% Output Argument:
%
%   Mdl - VAR(p) model with the following properties:
%
%       o Description - Model summary description
%       o SeriesNames - Response series names
%       o NumSeries   - Dimensionality of the response vector
%       o P           - Number of lagged responses of the VAR(p) model
%       o Constant    - Constant (intercept)
%       o AR          - Autoregressive coefficients
%       o Beta        - Regression coefficients
%       o Trend       - Linear time trend
%       o Covariance  - Innovations covariance matrix
%
% Notes:
%
%   o The short-hand syntax allows users to easily create model templates in
%     which the model dimensions are specified explicitly. The coefficients 
%     of these templates are then suitable for unrestricted parameter 
%     estimation, or may be updated via subsequent property assignment.
%
%     In contrast, the long-hand syntax is more powerful and allows users to
%     create models in which some or all of the coefficients are specified
%     and in which the model dimensions are inferred from the sizes of the
%     input coefficients.
%
%     In either syntax, once the model is created the number of time series 
%     (numseries) is fixed and cannot be changed.
%
%   o VAR models are most commonly expressed in the form of a difference
%     equation shown above. Specifying AR as a cell vector indicates that it 
%     is expressed in "difference equation format." In this convention, the 
%     sign and order of AR coefficients are determined exactly as encountered 
%     when reading the difference equation. Moreover, the difference equation 
%     shown above is expressed in reduced-form such that the coefficient of 
%     y(t) is assumed to be an identity matrix and is excluded from the cell
%     vector.
%
% See also ESTIMATE, FORECAST, INFER, SIMULATE.

% Copyright 2018 The MathWorks, Inc.

properties (GetAccess = public, SetAccess = private, Dependent)
  NumSeries                    % Dimensionality of the response vector (the number of responses)
  P                            % Number of lagged responses y(t) required to initialize the VAR model
end

properties (GetAccess = public, SetAccess = private, Dependent, Hidden)
  UnconditionalMean            % Unconditional (long-run) mean
end

properties (GetAccess = public, SetAccess = public, Dependent)
  Constant                     % Constant (intercept)
  AR                           % Auto-regressive coefficients of the P lagged responses
  Beta                         % Regression coefficients
  Trend                        % Linear time trend
  Covariance                   % Innovations covariance
  SeriesNames                  % K-element string vector of response series names
end

properties
  Description                  % Model description/summary string
end

properties (Access = private)
%
% The private properties are initialized to an empty state. Since the constructor 
% sets all private properties to some valid state in which K > 0, NON-EMPTY 
% private properties indicate that the model has already been created successfully.
%
% The only exception to this rule is the "FitInformation" property, which is non-empty 
% if the model has been estimated.
%
  PrivateConstant    = [];     % Constant (intercept)
  PrivateTrend       = [];     % Linear time trend
  PrivateBeta        = [];     % Regression coefficients
  PrivateCovariance  = [];     % Innovations covariance
  PrivateARLagOp     = {};     % Lag operator polynomial of the P autoregressive lagged responses
  PrivateSeriesNames = '';     % K-element string vector of response series names
  FitInformation     = [];     % Estimation information
end

methods    % Property GET & SET methods

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function K = get.NumSeries(Mdl)
     K = size(Mdl.PrivateConstant,1);
  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function P = get.P(Mdl)
     P = Mdl.PrivateARLagOp.Degree;
  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function Mdl = set.Description(Mdl, value)
     validateattributes(value, {'char'  'string'}, {'scalartext'}, '', 'Description');
     Mdl.Description = string(value);
  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function value = get.SeriesNames(Mdl)
     value = Mdl.PrivateSeriesNames;
  end

  function Mdl = set.SeriesNames(Mdl, value)
     K = Mdl.NumSeries;
     validateattributes(value, {'char' 'cell' 'string'}, {'vector'}, '', 'SeriesNames');
     if numel(value) ~= K
        error(message('econ:varm:varm:InconsistentSeriesSize', K))
     else
        Mdl.PrivateSeriesNames = string(value);
     end
  end
  
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function value = get.Constant(Mdl)
     value = Mdl.PrivateConstant;
  end
   
  function Mdl = set.Constant(Mdl, value)
     validateattributes(value, {'double'}, {'column'}, '', 'Constant');
     K = size(Mdl.PrivateConstant,1);
     if isscalar(value)
        Mdl.PrivateConstant = value(ones(K,1));
     else
        if size(value,1) ~= K
           error(message('econ:varm:varm:InconsistentConstantSize', K))
        end
        Mdl.PrivateConstant = value;
     end
     Mdl.FitInformation = [];  % Indicate that the model is NOT estimated.
  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function value = get.Trend(Mdl)
     value = Mdl.PrivateTrend;
  end
   
  function Mdl = set.Trend(Mdl, value)
     validateattributes(value, {'double'}, {'column'}, '', 'Trend');
     K                 = size(Mdl.PrivateTrend,1);
     updateDescription = Mdl.Description == getModelSummary(Mdl);
     if isscalar(value)
        Mdl.PrivateTrend = value(ones(K,1));
     else
        if size(value,1) ~= K
           error(message('econ:varm:varm:InconsistentTrendSize', K))
        end
        Mdl.PrivateTrend = value;
     end
     Mdl.FitInformation = [];  % Indicate that the model is NOT estimated.
     if updateDescription
        Mdl.Description = getModelSummary(Mdl);
     end
  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function value = get.Beta(Mdl)
     value = Mdl.PrivateBeta;
  end

  function Mdl = set.Beta(Mdl, value)
     validateattributes(value, {'double'}, {'2d'}, '', 'Beta');
     K                 = size(Mdl.PrivateBeta,1);
     updateDescription = Mdl.Description == getModelSummary(Mdl);
     if isscalar(value)
%
%       The following test allows uni-variate (i.e., 1-D) models to either
%       add a single predictor a model that has no predictors or remove 
%       predictors such that the resulting model has only a single predictor.
%       In the former case, a single predictor added to a 1-D model is NOT
%       a convenience interface designed to mimic scalar expansion, but
%       rather a direct, dimensionally-consistent assignment. 
%
%       The code below effects scalar expansion such that the Beta property 
%       is handled consistently as other properties (i.e., Constant, Trend). 
%       However, since the column dimension of Beta is determined by the 
%       number of predictors (i.e., it is NOT determined by the fundamental 
%       dimensions of the model), 1-D models are handled as a special case.
%
        if (K == 1) && (size(Mdl.PrivateBeta,2) <= 2)
           Mdl.PrivateBeta = value;
        else
           Mdl.PrivateBeta = value(ones(K,size(Mdl.PrivateBeta,2)));
        end
     else
        if size(value,1) ~= K
           error(message('econ:varm:varm:InconsistentBetaSize', K))
        end
        Mdl.PrivateBeta = value;
     end
     Mdl.FitInformation = [];  % Indicate that the model is NOT estimated.
     if updateDescription
        Mdl.Description = getModelSummary(Mdl);
     end
  end
  
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function value = get.Covariance(Mdl)
     value = Mdl.PrivateCovariance;
  end
   
  function Mdl = set.Covariance(Mdl, value)
     validateattributes(value, {'double'}, {'square'}, '', 'Covariance')
     K = size(Mdl.PrivateCovariance,1);
     if size(value,1) ~= K
        error(message('econ:varm:varm:InconsistentCovarianceSize', K, K))
     end
     temp              = value;
     temp(isnan(temp)) = 0;
     if ~issymmetric(temp)   % ensure symmetry
        error(message('econ:varm:varm:NonSymmetricCovariance', K, K))
     end
     if all(~isnan(value(:)))
        R = rank(value);
        if R ~= K            % ensure full-rank
           error(message('econ:varm:varm:NonFullRankCovariance', R, K))
        end
     end
     Mdl.PrivateCovariance = value;
     Mdl.FitInformation    = [];  % Indicate that the model is NOT estimated.
  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function coefficients = get.AR(Mdl)
     if Mdl.PrivateARLagOp.Degree > 0
        coefficients = toCellArray(reflect(Mdl.PrivateARLagOp));
        coefficients = coefficients(2:end);   % Coefficients at lags 1, 2, ...
     else
        coefficients = {};
     end
  end

  function Mdl = set.AR(Mdl, value)
     if isempty(value)
        validateattributes(value, {'cell'}, {}, '', 'AR');
     else
        validateattributes(value, {'cell'}, {'vector'}, '', 'AR');
     end
     K                 = size(Mdl.PrivateConstant,1);
     updateDescription = Mdl.Description == getModelSummary(Mdl);
     for i = 1:numel(value)
         if isempty(value{i})
            value{i} = zeros(K);
         elseif isscalar(value{i})
            value{i} = value{i}(ones(K));
         else
            if (size(value{i},1) ~= K) || (size(value{i},2) ~= K)
               error(message('econ:varm:varm:InconsistentARSize', K, K))
            end
         end
     end
     Mdl.PrivateARLagOp = LagOp([Mdl.PrivateARLagOp.Coefficients{0} cellfun(@uminus, value(:)', 'uniformoutput', false)]);
     Mdl.FitInformation = [];  % Indicate that the model is NOT estimated.
     if updateDescription
        Mdl.Description = getModelSummary(Mdl);
     end
  end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

  function average = get.UnconditionalMean(Mdl)
     AR              = toCellArray(Mdl.PrivateARLagOp);                         % Get AR coefficients at all lags
     nPredictors     = size(Mdl.PrivateBeta,2);                                 % # of predictors
     isTrendIncluded = any(Mdl.PrivateTrend) || any(isnan(Mdl.PrivateTrend));   % Is a linear time trend included?
     if isStable(Mdl.PrivateARLagOp) && ~isTrendIncluded && (nPredictors == 0)     
        summation = AR{1};
        for Lag = 2:numel(AR)
            summation = summation + AR{Lag};              % A(L = 1) = I + A1 + A2 + ... + Ap)
        end
        try
           average = summation \ Mdl.PrivateConstant;
        catch
           average = nan(size(Mdl.PrivateConstant,1),1);
        end
     else
        average = nan(size(Mdl.PrivateConstant,1),1);
     end
  end

end  % METHODS Block


methods (Access = public)

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function Mdl = varm(varargin)
%VARM Construct a VAR(p) model

if nargin == 0              % Is it the trivial syntax Mdl = varm()?
% 
%  MATLAB classes should construct scalar objects in default states in the 
%  absence of input arguments. In this case, the default constructor syntax 
%  creates a 1-D VAR(0) model (i.e., a simple, univariate random walk) with 
%  an undefined constant and covariance (i.e., NaN).
%
   K                     = 1;
   Mdl.PrivateConstant   = nan(K,1);
   Mdl.PrivateBeta       = nan(K,0);
   Mdl.PrivateTrend      = zeros(K,1);
   Mdl.PrivateCovariance = nan(K);
   Mdl.PrivateARLagOp    = LagOp(eye(K));
   Mdl.Description       = getModelSummary(Mdl);
   Mdl.SeriesNames       = "Y";
   return
end

%
% Validate input parameters.
%

if isnumeric(varargin{1}) && (nargin == 2)  % Is it the short-hand Mdl = varm(K,P) syntax?

%
%  Validate short-hand syntax.
%
   parser = inputParser;
   parser.addRequired('numseries', @(x) validateattributes(x, {'double'}, {'scalar' 'positive'    'integer'}, '', 'numseries'));
   parser.addRequired('numlags'  , @(x) validateattributes(x, {'double'}, {'scalar' 'nonnegative' 'integer'}, '', 'numlags'));
   parser.parse(varargin{:});

   K = parser.Results.numseries;
   P = parser.Results.numlags ;
%
%  Initialize model properties. 
%
   Mdl.PrivateARLagOp    = LagOp([{eye(K)} repmat({nan(K)},1,P)]);
   Mdl.PrivateConstant   = nan(K,1);
   Mdl.PrivateTrend      = zeros(K,1);
   Mdl.PrivateBeta       = nan(K,0);
   Mdl.PrivateCovariance = nan(K,K);
   Mdl.Description       = getModelSummary(Mdl);
   Mdl.SeriesNames       = repmat("Y", 1, K);

   for i = 1:K
       Mdl.SeriesNames(i) = Mdl.SeriesNames(i) + i;
   end

elseif ischar(varargin{1}) || ( isstring(varargin{1}) && isscalar(varargin{1}) )  % Is it the long-hand syntax of N-V pairs?

%
%  Validate long-hand syntax of parameter name-value pairs.
%

   try
     Mdl = validate(Mdl, varargin{:});
   catch exception
     exception.throwAsCaller();
   end

else
   error(message('econ:varm:varm:UnsupportedConstructorSyntax'))
end

end % Constructor

end % Methods (Access = public)


methods (Hidden)

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function polynomial = getLagOp(Mdl, name)
%GETLAGOP Get VAR model autoregressive lag operator polynomial
%
% Syntax:
%
%   polynomial = getLagOp(Mdl,name)
%
% Description:
%
%   Return the underlying lag operator polynomial (LagOp) of the autoregressive 
%   coefficients associated with lagged responses. The degree of the autoregressive 
%   polynomial is P.
%
% Input Arguments:
%
%   Mdl - VAR(p) model whose autoregressive polynomial is requested.
%
%   name - The name of the polynomial (case-insensitive character string).
%     The available option is 'AR'.
%
% Output Arguments:
%
%   polynomial - Autoregressive lag operator polynomial of the VAR(p) model.

   switch upper(name)                    % Name of the component polynomial

      case upper('AR')                   % Autoregressive coefficient polynomial
         polynomial = Mdl.PrivateARLagOp;

      otherwise
         error(message('econ:varm:varm:InvalidPolynomial'))
   end

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function varargout = getModelSummary(Mdl)
%
% Create the model summary string used by the "displayScalarObject" 
% and "summarize" methods.
%

K  = Mdl.NumSeries;
ar = [Mdl.AR{:}];

%
% Create & print the model summary string.
%

if any(isnan(ar(:))) || (Mdl.P == 0)
   summary = '   ';
else
   if isStable(Mdl.PrivateARLagOp)
      summary = '   AR-Stationary ';
   else
      summary = '   AR-Nonstationary ';
   end
end

if isempty(Mdl.PrivateBeta) || all(Mdl.PrivateBeta(:) == 0)
   summary = [summary  '%d-Dimensional VAR(%d) Model'];
else
   if size(Mdl.PrivateBeta,2) == 1
      summary = [summary  '%d-Dimensional VARX(%d) Model with %d Predictor'];
   else
      summary = [summary  '%d-Dimensional VARX(%d) Model with %d Predictors'];
   end
end

if any(Mdl.PrivateTrend) || any(isnan(Mdl.PrivateTrend))
   if isempty(Mdl.PrivateBeta)
      summary = [summary  ' with Linear Time Trend'];
   else
      summary = [summary  ' and Linear Time Trend'];
   end
end

if nargout > 0
   if isempty(Mdl.PrivateBeta) || all(Mdl.PrivateBeta(:) == 0)
      varargout = {string(strtrim(sprintf(summary, K, Mdl.P)))};
   else
      varargout = {string(strtrim(sprintf(summary, K, Mdl.P, size(Mdl.PrivateBeta,2))))};
   end 
else
   if isempty(Mdl.PrivateBeta) || all(Mdl.PrivateBeta(:) == 0)
      fprintf(['<strong>' summary '</strong>\n'], K, Mdl.P);
   else
      fprintf(['<strong>' summary '</strong>\n'], K, Mdl.P, size(Mdl.PrivateBeta,2));
   end 
end

end % Get Model Summary

end % Methods (Hidden)


methods (Access = protected)

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function Mdl = validate(Mdl, varargin)
%
% Validate a variable-length list (long-hand syntax) of parameter name-value 
% pairs passed to the constructor.
%

parser = inputParser;
parser.addParameter('Constant'   , Mdl.PrivateConstant  , @(x) validateattributes(x, {'double'}       , {'column'}                     , '', 'Constant'));
parser.addParameter('AR'         , Mdl.PrivateARLagOp   , @(x) validateattributes(x, {'cell'}         , {}                             , '', 'AR'));
parser.addParameter('Lags'       , []                   , @(x) validateattributes(x, {'double'}       , {'vector' 'integer' 'positive'}, '', 'Autoregressive Lags'));
parser.addParameter('Trend'      , Mdl.PrivateTrend     , @(x) validateattributes(x, {'double'}       , {'column'}                     , '', 'Time Trend'));
parser.addParameter('Beta'       , Mdl.PrivateBeta      , @(x) validateattributes(x, {'double'}       , {'2d'}                         , '', 'Beta'));
parser.addParameter('Covariance' , Mdl.PrivateCovariance, @(x) validateattributes(x, {'double'}       , {'square'}                     , '', 'Covariance'));
parser.addParameter('Description', ''                   , @(x) validateattributes(x, {'char' 'string'}, {'scalartext'}                 , '', 'Description'));
parser.addParameter('SeriesNames', ''                   , @(x) validateattributes(x, {'char' 'cell' 'string'}, {'vector'}              , '', 'SeriesNames'));

parser.parse(varargin{:});
   
constant    = parser.Results.Constant;
ar          = parser.Results.AR;
lags        = parser.Results.Lags;
trend       = parser.Results.Trend;
beta        = parser.Results.Beta;
covariance  = parser.Results.Covariance;
description = string(parser.Results.Description);
series      = string(parser.Results.SeriesNames);

if ~(isempty(ar) || isvector(ar))
   error(message('econ:varm:varm:InvalidAutoregressiveVector'))
end

%
% Initialize the number of series (K = numseries) to a "has yet to be determined" state.
%

K = NaN;

%
% Create a series of flags and determine the number of series (K = numseries), 
% which is a fundamental dimension of VAR models and must be determined to 
% create a valid model. 
%
% In what follows, K is determined by the first parameter checked whose 
% dimension is associated with K.
%

isConstantSpecified   = ~any(strcmpi('Constant'  , parser.UsingDefaults));
isARSpecified         = ~any(strcmpi('AR'        , parser.UsingDefaults));
isTrendSpecified      = ~any(strcmpi('Trend'     , parser.UsingDefaults));
isBetaSpecified       = ~any(strcmpi('Beta'      , parser.UsingDefaults));
isCovarianceSpecified = ~any(strcmpi('Covariance', parser.UsingDefaults));

if isConstantSpecified
   K = size(constant,1);   
elseif isARSpecified && (numel(ar) > 0) 
   K = size(ar{1},1);   
elseif isTrendSpecified 
   K = size(trend,1);   
elseif isBetaSpecified 
   K = size(beta,1);   
elseif isCovarianceSpecified
   K = size(covariance,1); 
end
   
if isnan(K)
   error(message('econ:varm:varm:UndeterminedSeriesSize'))
end

if K == 0
   error(message('econ:varm:varm:InsufficientSeriesSize'))
end

%
% Now that the number of series (K = numseries) is known, enforce dimensional 
% consistency between the various model parameters and assign defaults if necessary.
%

if isConstantSpecified
   n = size(constant,1);
   if n ~= K
      error(message('econ:varm:varm:InconsistentConstant', n, K))
   end
else
   constant = nan(K,1);
end

if isTrendSpecified
   n = size(trend,1); 
   if n ~= K
      error(message('econ:varm:varm:InconsistentTrend', n, K))
   end
else
   trend = zeros(K,1);
end

if isBetaSpecified
   n = size(beta,1); 
   if n ~= K
      error(message('econ:varm:varm:InconsistentBeta', n, size(beta,2), K))
   end
else
   beta = nan(K,0);
end

if isCovarianceSpecified
   [n,m]= size(covariance);
   if (n ~= K) && (m ~= K)
      error(message('econ:varm:varm:InconsistentCovariance', n, m, K, K))
   end
   temp              = covariance;
   temp(isnan(temp)) = 0;
   if ~issymmetric(temp)             % ensure symmetry
      error(message('econ:varm:varm:NonSymmetricCovariance', K, K))
   end
   if all(~isnan(covariance(:)))
      R = rank(covariance); 
      if R ~= K                      % ensure full-rank
         error(message('econ:varm:varm:InconsistentCovarianceRank', R, K, K))
      end
   end
else
   covariance = nan(K,K);
end
  
%
% Check autoregressive coefficients and corresponding lags.
%
% In the following, if the user specified coefficients without lags, then 
% derive lags from the coefficients. Similarly, if the user specified lags 
% without coefficients, then assign NaNs to coefficients consistent with 
% the number of lags.
%

ar = ar(:)';

for i = 1:numel(ar)
    if (size(ar{i},1) ~= K) || (size(ar{i},2) ~= K)
       error(message('econ:varm:varm:InconsistentARDimensions', K, K))
    end
end

if any(strcmpi('Lags', parser.UsingDefaults))
   lags = 0:numel(ar);
else
   if any(strcmpi('AR', parser.UsingDefaults))
      ar = repmat({nan(K)}, 1, numel(lags));
   end

   if numel(lags) ~= numel(ar)
      error(message('econ:varm:varm:ARLagInconsistency'))
   end

   uniqueLags = unique(lags, 'first');

   if any(lags == 0) || (numel(lags) ~= numel(uniqueLags))
      error(message('econ:varm:varm:InconsistentARLags'))
   end
   lags = [0 lags];
end

%
% Create the lag operator polynomial associated with the underlying AR
% coefficients and initialize the remaining properties. Note that the AR 
% coefficients are specified as a cell array and therefore indicate a 
% difference equation expressed in reduced-form (i.e., the zero-lag AR
% coefficient is guaranteed to be I = eye(K)).
%

if isempty(Mdl.PrivateARLagOp)
   Mdl.PrivateARLagOp = LagOp([eye(K) cellfun(@uminus, ar, 'uniformoutput', false)], 'Lags', lags);
end

Mdl.PrivateConstant   = constant;
Mdl.PrivateTrend      = trend;
Mdl.PrivateBeta       = beta;
Mdl.PrivateCovariance = covariance;

%
% Assign or auto-generate defaults for the model description and response
% series names.
%

if ~any(strcmpi('Description', parser.UsingDefaults))
   Mdl.Description = description;
else
   Mdl.Description = getModelSummary(Mdl);
end

if ~any(strcmpi('SeriesNames', parser.UsingDefaults))
   if numel(series) == K
      Mdl.SeriesNames = series;
   else
      error(message('econ:varm:varm:MismatchedSeriesNames'))
   end
else
   Mdl.SeriesNames = repmat("Y", 1, K);
   for i = 1:K
       Mdl.SeriesNames(i) = Mdl.SeriesNames(i) + i;
   end
end

end  % Validate

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function Mdl = setLagOp(Mdl, name, polynomial)
%SETLAGOP Set VAR model autoregressive lag operator polynomial
%
% Syntax:
%
%   Mdl = setLagOp(Mdl,name,polynomial)
%
% Description:
%
%   Set the underlying lag operator polynomial (LagOp) of the autoregressive 
%   coefficients associated with lagged responses.
%
% Input Arguments:
%
%   Mdl - VAR model whose component polynomial is updated.
%
%   name - The name of the polynomial (case-insensitive character string).
%     The available option is 'AR'.
%
%   polynomial - Updated autoregressive lag operator polynomial.
%
% Output Arguments:
%
%   Mdl - VAR(p) model whose autoregressive polynomial is updated.

   switch upper(name)                        % Name of the component polynomial

      case upper('AR')                       % Autoregressive coefficient polynomial
         if Mdl.PrivateARLagOp.Dimension ~= polynomial.Dimension
            error(message('econ:varm:varm:InvalidPolynomialDimension', Mdl.NumSeries))
         end
         Mdl.PrivateARLagOp = polynomial;

      otherwise
         error(message('econ:varm:varm:InvalidPolynomial'))
   end

end

%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%

function displayScalarObject(Mdl)

%
% Override matlab.mixin.CustomDisplay.displayScalarObject method.
%

maxToPrint = 3;  % Maximum # of coefficients to explicitly print
K          = Mdl.NumSeries;
x          = char(215);

disp(matlab.mixin.CustomDisplay.getSimpleHeader(Mdl));

%
% Print model description and response series names.
%

spaces = '                    ';
fprintf([spaces(1:5) 'Description: "%s"\n'], Mdl.Description)

if K <= maxToPrint
   format = repmat(' "%s" ', 1, K);
   fprintf(string(spaces(1:5)) + 'SeriesNames: ' + format(2:end) + '\n', Mdl.SeriesNames(1:K));
else
   format = repmat(' "%s" ', 1, maxToPrint);
   fprintf(string(spaces(1:5)) + 'SeriesNames: ' + format(2:end) + ' ... and %d more\n', Mdl.SeriesNames(1:maxToPrint), K - maxToPrint);
end

fprintf([spaces(1:7) 'NumSeries: %d\n'], Mdl.NumSeries)
fprintf([spaces(1:7) '        P: %d\n'], Mdl.P)

%
% Print the K-by-1 constant.
%

if K == 1
   fprintf([spaces(1:8) 'Constant: %g\n'], Mdl.PrivateConstant)
else
   if all(Mdl.PrivateConstant == 0)
      fprintf([spaces(1:8) 'Constant: [%d' x '%d vector of zeros]\n'], K, 1)
   elseif all(isnan(Mdl.PrivateConstant))
      fprintf([spaces(1:8) 'Constant: [%d' x '%d vector of NaNs]\n'], K, 1)
   else
      if K <= maxToPrint
         format = repmat(' %g', 1, K);
         fprintf([spaces(1:8) 'Constant: [' format(2:end) ']''\n'], Mdl.PrivateConstant(1:K)')
      else
         format = repmat(' %g', 1, maxToPrint);
         fprintf([spaces(1:8) 'Constant: [' format(2:end) ' ... and %d more]''\n'], Mdl.PrivateConstant(1:maxToPrint)', K - maxToPrint)
      end
   end
end

% 
% Print K-by-K AR coefficients.
%

L  = Mdl.PrivateARLagOp.Lags;
L  = L(L > 0);
N  = numel(L);                   % # of non-zero lags included 
ar = [Mdl.AR{:}];

if N == 0
   fprintf([spaces(1:14) 'AR: {}\n'])
elseif N <= maxToPrint
   format = repmat(' %g', 1, N);
   if K == 1
      if N == 1
         fprintf([spaces(1:14) 'AR: {' format(2:end) '} at lag [' format(2:end) ']\n'], ar(L), L)
      else
         fprintf([spaces(1:14) 'AR: {' format(2:end) '} at lags [' format(2:end) ']\n'], ar(L), L)
      end
   else
      if N == 1
         if all(isnan(ar(:)))
            fprintf([spaces(1:14) 'AR: {%d' x '%d matrix of NaNs} at lag [' format(2:end) ']\n'], K, K, L)
         else
            fprintf([spaces(1:14) 'AR: {%d' x '%d matrix} at lag [' format(2:end) ']\n'], K, K, L)
         end
      else
         if all(isnan(ar(:)))
            fprintf([spaces(1:14) 'AR: {%d' x '%d matrices of NaNs} at lags [' format(2:end) ']\n'], K, K, L)
         else
            fprintf([spaces(1:14) 'AR: {%d' x '%d matrices} at lags [' format(2:end) ']\n'], K, K, L)
         end
      end
   end
else
   format = repmat(' %g', 1, maxToPrint);
   if K == 1
      if all(isnan(ar(:)))
         fprintf([spaces(1:14) 'AR: {NaNs} at lags [' format(2:end) ' ... and %d more]\n'], L(1:maxToPrint), N - maxToPrint)
      else
         fprintf([spaces(1:14) 'AR: {' format(2:end) ' ... and %d more} at lags [' format(2:end) ' ... and %d more]\n'], ar(L(1:maxToPrint)), N - maxToPrint, L(1:maxToPrint), N - maxToPrint)
      end
   else
      if all(isnan(ar(:)))
         fprintf([spaces(1:14) 'AR: {%d' x '%d matrices of NaNs} at lags [' format(2:end) ' ... and %d more]\n'], K, K, L(1:maxToPrint), N - maxToPrint)
      else
         fprintf([spaces(1:14) 'AR: {%d' x '%d matrices} at lags [' format(2:end) ' ... and %d more]\n'], K, K, L(1:maxToPrint), N - maxToPrint)
      end
   end
end

%
% Print the K-by-1 trend.
%

if K == 1
   fprintf([spaces(1:11) 'Trend: %g\n'], Mdl.PrivateTrend)
else
   if all(Mdl.PrivateTrend == 0)
      fprintf([spaces(1:11) 'Trend: [%d' x '%d vector of zeros]\n'], K, 1)
   elseif all(isnan(Mdl.PrivateTrend))
      fprintf([spaces(1:11) 'Trend: [%d' x '%d vector of NaNs]\n'], K, 1)
   else
      if K <= maxToPrint
         format = repmat(' %g', 1, K);
         fprintf([spaces(1:11) 'Trend: [' format(2:end) ']''\n'], Mdl.PrivateTrend(1:K)')
      else
         format = repmat(' %g', 1, maxToPrint);
         fprintf([spaces(1:11) 'Trend: [' format(2:end) ' ... and %d more]''\n'], Mdl.PrivateTrend(1:maxToPrint)', K - maxToPrint)
      end
   end
end

%
% Print the K-by-N regression coefficients.
%

N = size(Mdl.PrivateBeta,2);    % # of predictors

if K == 1
   if N <= maxToPrint
      format = repmat(' %g', 1, N);
      if N == 0
        fprintf([spaces(1:12) 'Beta: [%d' x '%d matrix]\n'], K, N)
      else
        fprintf([spaces(1:12) 'Beta: [' format(2:end) ']\n'], Mdl.PrivateBeta(1:N)')
      end
   else
      if all(Mdl.PrivateBeta(:) == 0)
         fprintf([spaces(1:12) 'Beta: [%d' x '%d matrix of zeros]\n'], K, N)
      elseif all(isnan(Mdl.PrivateBeta(:)))
         fprintf([spaces(1:12) 'Beta: [%d' x '%d matrix of NaNs]\n'], K, N)
      else
         format = repmat(' %g', 1, maxToPrint);
         fprintf([spaces(1:12) 'Beta: [' format(2:end) ' ... and %d more]\n'], Mdl.PrivateBeta(1:maxToPrint)', N - maxToPrint)
      end
   end      
else
   if isempty(Mdl.PrivateBeta)
      fprintf([spaces(1:12) 'Beta: [%d' x '%d matrix]\n'], K, N)
   elseif all(Mdl.PrivateBeta(:) == 0)
      fprintf([spaces(1:12) 'Beta: [%d' x '%d matrix of zeros]\n'], K, N)
   elseif all(isnan(Mdl.PrivateBeta(:)))
      fprintf([spaces(1:12) 'Beta: [%d' x '%d matrix of NaNs]\n'], K, N)
   else
      fprintf([spaces(1:12) 'Beta: [%d' x '%d matrix]\n'], K, N)
   end
end

%
% Print the K-by-K covariance.
%

if K == 1
   fprintf([spaces(1:6) 'Covariance: %g\n'], Mdl.PrivateCovariance)
else
   if isdiag(Mdl.PrivateCovariance)
      fprintf([spaces(1:6) 'Covariance: [%d' x '%d diagonal matrix]\n'], K, K)
   else
      if all(isnan(Mdl.PrivateCovariance(:)))
         fprintf([spaces(1:6) 'Covariance: [%d' x '%d matrix of NaNs]\n'], K, K)
      else
         fprintf([spaces(1:6) 'Covariance: [%d' x '%d matrix]\n'], K, K)
      end
   end
end

%
% Print the model summary string if needed.
%

summary = getModelSummary(Mdl);

if Mdl.Description ~= summary
   disp(' ')
   fprintf("   " + '<strong>' + summary + '</strong>\n');
end

end % Display Scalar Object

end % Methods (Protected)

end % Class definition