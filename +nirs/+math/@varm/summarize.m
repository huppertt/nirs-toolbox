function varargout = summarize(Mdl)
%SUMMARIZE Summarize vector autoregression (VAR) model estimation results
%
% Syntax:
%
%   summarize(Mdl)
%   results = summarize(Mdl)
%
% Description:
%
%   Given an estimated vector autoregression (VAR) model, summarize the 
%   estimation results, including a table of model parameter values, standard
%   errors, t statistics, and p-values. The summary also includes various
%   statistics such as loglikelihood, AIC/BIC, and innovations correlation 
%   matrix.
%
%   When called without an output (first syntax), summary information is 
%   printed to the MATLAB command window. When called with an output, no 
%   information is printed, but the summary is returned in the results
%   structure.
%
%   If the model is not estimated, then it is simply displayed to the command
%   window when called without an output; otherwise the input model is 
%   returned.
%
% Input Arguments:
%
%   Mdl - VAR model created by the VARM constructor or VARM/ESTIMATE method.
%
% Output Arguments:
%
%   results - For estimated models, a structure of estimation summary 
%     information with the following fields:
%
%     o Description            - Model summary description
%     o SampleSize             - Effective sample size
%     o NumEstimatedParameters - Number of estimated parameters
%     o LogLikelihood          - Loglikelihood value
%     o AIC                    - Akaike information criterion
%     o BIC                    - Bayesian information criterion
%     o Table                  - Table of parameter values, standard errors, t 
%                                statistics (value divided by standard error),
%                                and p-values (assuming normality)
%     o Covariance             - Residual covariance matrix (MLE)
%     o Correlation            - Residual correlation matrix
%
%     If the model is not estimated, then the input model (Mdl) is returned.
%
% See also ESTIMATE.

% Copyright 2017 The MathWorks, Inc.

if numel(Mdl) > 1
   error(message('econ:varm:summarize:NonScalarModel'))
end

%
% If the model has not been estimated, then simply call the display method
% or return the input model.
%

if isempty(Mdl.FitInformation)
   if nargout > 0
      varargout = {Mdl};
   else
      displayScalarObject(Mdl)
   end
   return
end

%
% Initialize some basic model information.
%

K                  = Mdl.NumSeries;
P                  = Mdl.P;
Lags               = Mdl.PrivateARLagOp.Lags;  % Lags included in the AR polynomial (includes zero)
nAR                = sum(Lags > 0); 
Lags               = Lags(Lags > 0);    
isConstantIncluded = any(Mdl.PrivateConstant);
isTrendIncluded    = any(Mdl.PrivateTrend);
nPredictors        = size(Mdl.PrivateBeta,2);

%
% Compute the actual number of parameters included in the model, which determines
% number of rows in the parameter Table.
%

nParameters  = (K * isConstantIncluded) + (nAR * K * K) + (nPredictors * K) + (K * isTrendIncluded);

%
% Create the row names and column information of the summary table.
%

rowNames      = cell(nParameters,1);
Value         = zeros(nParameters,1);
StandardError = zeros(nParameters,1);
iRow          = 1;

if isConstantIncluded      % Constant
   for i = 1:K
       rowNames{iRow}      = sprintf('Constant(%d)', i);
       Value(iRow)         = Mdl.PrivateConstant(i);
       StandardError(iRow) = Mdl.FitInformation.StandardErrors.Constant(i);
       iRow                = iRow + 1;
   end
end

for i = 1:P                % AR terms
    if any(i == Lags)
       for i2 = 1:K
           for i1 = 1:K
               rowNames{iRow}      = sprintf('AR{%d}(%d,%d)', i, i1, i2);
               Value(iRow)         = Mdl.AR{i}(i1,i2);
               StandardError(iRow) = Mdl.FitInformation.StandardErrors.AR{i}(i1,i2);
               iRow                = iRow + 1;
           end
        end
    end
end

for i2 = 1:nPredictors      % Predictors
    for i1 = 1:K
        rowNames{iRow}      = sprintf('Beta(%d,%d)', i1, i2);
        Value(iRow)         = Mdl.PrivateBeta(i1,i2);
        StandardError(iRow) = Mdl.FitInformation.StandardErrors.Beta(i1,i2);
        iRow                = iRow + 1;
    end
end

if isTrendIncluded         % Linear time rend
   for i = 1:K
       rowNames{iRow}      = sprintf('Trend(%d)', i);
       Value(iRow)         = Mdl.PrivateTrend(i);
       StandardError(iRow) = Mdl.FitInformation.StandardErrors.Trend(i);
       iRow                = iRow + 1;
   end
end

%
% Compute corresponding t-statistics & p-values and create the summary table.
%

tStatistic = Value ./ StandardError;
pValue     = 2 * (normcdf(-abs(tStatistic)));
Table      = table(Value, StandardError, tStatistic, pValue, 'RowNames', rowNames);
Table.Properties.VariableNames = {'Value' 'StandardError' 'TStatistic' 'PValue'};

%
% Compute additional statistics.
%

if Mdl.FitInformation.NumEstimatedParameters > 0
  [AIC,BIC] = aicbic(Mdl.FitInformation.LogLikelihood, Mdl.FitInformation.NumEstimatedParameters, Mdl.FitInformation.SampleSize);
else
   AIC = -2 * Mdl.FitInformation.LogLikelihood;
   BIC = -2 * Mdl.FitInformation.LogLikelihood;
end

%
% Compute innovations correlation matrix (see COV2CORR in Financial Toolbox).
%

sigma                                 = sqrt(diag(Mdl.PrivateCovariance));
correlation                           = Mdl.PrivateCovariance ./ (sigma * sigma');
correlation(sub2ind([K K], 1:K, 1:K)) = 1; % Force exact ones along the main diagonal.

%
% Now print or return the summary information.
%

if nargout > 0
   output.Description            = Mdl.Description;
   output.SampleSize             = Mdl.FitInformation.SampleSize;
   output.NumEstimatedParameters = Mdl.FitInformation.NumEstimatedParameters;
   output.LogLikelihood          = Mdl.FitInformation.LogLikelihood;
   output.AIC                    = AIC;
   output.BIC                    = BIC;
   output.Table                  = Table;
   output.Covariance             = Mdl.PrivateCovariance;
   output.Correlation            = correlation;
   varargout                     = {output};
else
   disp(' ')
   fprintf("   " + '<strong>' + Mdl.Description + '</strong>\n');
   disp(' ')
   fprintf('    Effective Sample Size: %d\n', Mdl.FitInformation.SampleSize)
   fprintf('    Number of Estimated Parameters: %d\n', Mdl.FitInformation.NumEstimatedParameters)
   fprintf('    LogLikelihood: %g\n', Mdl.FitInformation.LogLikelihood)
   fprintf('    AIC: %g\n', AIC)
   fprintf('    BIC: %g\n', BIC)
   disp(' ')
   disp(Table)
   disp(' ')
   fprintf('<strong>   Innovations Covariance Matrix:</strong>\n')
   disp(Mdl.PrivateCovariance)
   disp(' ')
   fprintf('<strong>   Innovations Correlation Matrix:</strong>\n')
   disp(correlation)
end

end               