function [p,stat] = dwtest(model,option,tail)
%DWTEST Durbin-Watson test for autocorrelation.
%    [P,DW]=DWTEST(LM) performs a Durbin-Watson test on the residuals from
%    the LinearModel LM.  P is the computed p-value for the test, and DW is
%    the Durbin-Watson statistic.  The Durbin-Watson test is used to test
%    if the residuals are uncorrelated, against the alternative that there
%    is autocorrelation among them.
%
%    [...]=DWTEST(LM,METHOD) specifies the method to be used in
%    computing the p-value.  METHOD can be either of the following:
%       'exact'        Calculate an exact p-value using the PAN algorithm
%                      (default if the sample size is less than 400).
%       'approximate'  Calculate the p-value using a normal approximation
%                      (default if the sample size is 400 or larger).
%
%    [...]=DWTEST(LM,METHOD,TAIL) performs the test against the
%    alternative hypothesis specified by TAIL:
%       'both'   "serial correlation is not 0" (two-tailed test, default)
%       'right'  "serial correlation is greater than 0" (right-tailed test)
%       'left'   "serial correlation is less than 0" (left-tailed test)
%
%    Example:
%       % There is significant autocorrelation in the residuals from a
%       % straight-line fit to the census data
%       load census
%       lm = fitlm(cdate,pop);
%       p = dwtest(lm)
%
%       % Residuals from a quadratic fit have reduced autcorrelation, but
%       % it is still significant
%       lm = fitlm(cdate,pop,'quadratic');
%       p = dwtest(lm)
%
%   See also LinearModel, Residuals.

%   Reference:
%   J. Durbin & G.S. Watson (1950), Testing for Serial Correlation in Least
%   Squares Regression I. Biometrika (37), 409-428.
%
%   R.W. Farebrother (1980), Pan's Procedure for the Tail Probabilities of
%   the Durbin-Watson Statistic. Applied Statistics, 29, 224-227.
%
%   Copyright 2011-2015 The MathWorks, Inc.


if nargin < 2
    if model.NumObservations < 400
        option = 'exact';
    else
        option = 'approximate';
    end;
end;
if nargin < 3
    tail = 'both'; % default test is two sided
end;

subset = model.ObservationInfo.Subset;
r = model.Residuals.Raw(subset);
stat = sum(diff(r).^2)/sum(r.^2); % durbin-watson statistic

% This calls the function of Pan algorithm/normal approximation
% Recall that the distribution of dw depends on the design matrix
% in the regression.
pdw = dfswitchyard('pvaluedw',stat,model.design_r,option);

% p-value depends on the alternative hypothesis.
switch lower(tail)
    case 'both'
        p = 2*min(pdw, 1-pdw);
    case 'left'
        p = 1-pdw; % *** compute upper tail directly?
    case 'right'
        p = pdw;
end

