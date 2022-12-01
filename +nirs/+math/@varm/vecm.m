function VEC = vecm(VAR)
%VECM Vector autoregression (VAR) model to vector error correction (VEC) model
%
% Syntax:
%
%   VEC = vecm(VAR)
%
% Description:
%
%   Given a vector autoregression (VAR) model, convert to an equivalent
%   vector error correction (VEC) model. A VAR(p) model is given by the 
%   difference equation
%
%   y(t) = c + B1*y(t-1) + ... + Bp*y(t-p) + D*x(t) + T*t + e(t)
%
%   for responses y(t), predictors x(t), and innovations e(t).
%
%   Using lag operator notation such that Ly(t) = y(t-1), the equivalent
%   VEC(p-1) model with q = p-1 lagged changes (1-L)y(t) is given by the
%   difference equation
%
%   (1-L)y(t) = c + C*y(t-1) + A1*(1-L)y(t-1) + ... + Aq*(1-L)y(t-q)
%                 + D*x(t)   + T*t + e(t)
%
%   where C is the impact (error correction) coefficient.
%
% Input:
%
%   VAR - Vector autoregression model, created by the VARM constructor or
%      the VARM/ESTIMATE method.
%
% Output:
%
%   VEC - Vector error correction model equivalent to the input VAR model.
%
% References:
%
%   [1] Johansen, S. Likelihood-Based Inference in Cointegrated Vector
%       Autoregressive Models. Oxford: Oxford University Press, 1995.
%
%   [2] Juselius, K. The Cointegrated VAR Model. Oxford: Oxford University 
%       Press, 2006.
%
%   [3] Hamilton, J. D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994.
%
%   [4] Lutkepophl, H. New Introduction to Multiple Time Series Analysis.
%       Springer-Verlag, 2007.
%
% See also VAR2VEC.

% Copyright 2017 The MathWorks, Inc.

if numel(VAR) > 1
   error(message('econ:varm:vecm:NonScalarModel'))
end

%
% Construct a VAR model & assign properties.
%

[shortRun,impact] = var2vec(VAR.AR);

VEC = vecm('Constant'   , VAR.PrivateConstant, 'Impact'     , impact               , ...
           'ShortRun'   , shortRun           , 'Beta'       , VAR.PrivateBeta      , ...
           'Trend'      , VAR.PrivateTrend   , 'Covariance' , VAR.PrivateCovariance, ...
           'Description', VAR.Description    , 'SeriesNames', VAR.SeriesNames      );
%
% Update the output VEC model description only if the input VAR model
% description was auto-generated. Otherwise, retain the description specified 
% by the user.
%

if VAR.Description == getModelSummary(VAR)
   VEC.Description = getModelSummary(VEC);
end

end
