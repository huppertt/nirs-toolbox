function t = tstats(coefs,se,nobs,coefnames)
%TSTATS Table of t statistics for a regression model fit.
%   TBL = TSTATS(COEFS,SE,NOBS) returns a table of t statistics for the
%   estimated coefficients in a regression fit (typically a linear or
%   nonlinear regression).  TBL is a dataset array containing entries (rows)
%   for each coefficient in the fitted model, and variables (columns) for the
%   coefficient estimates, their standard errors, the t statistics, and the
%   2-sided p-values.
%
%   COEFS is a vector containing the estimated regression coefficients, SE is
%   a vector containing their standard errors, and NOBS is the number of
%   observations used in the fit.
%
%   TBL = TSTATS(COEFS,SE,NOBS,COEFNAMES) uses COEFNAMES as the names for each
%   coefficient in the fit.
%
%   The p-values are based on the assumption that the estimated coefficients
%   in the regression are approximately unbiased and normally distributed, and
%   that the squares of the standard errors are approximately chi-squared
%   distributed.

%   Copyright 2011-2014 The MathWorks, Inc.

if ~internal.stats.isScalarInt(nobs,1)
    error(message('stats:classreg:regr:modelutils:BadNumObs'));
end

coefs = coefs(:);
ncoefs = numel(coefs);
se = se(:);
if numel(se) ~= ncoefs
    error(message('stats:classreg:regr:modelutils:CoefSEDifferentLength'));
end

if nargin < 4
    coefnames = internal.stats.numberedNames('b',1:ncoefs);
elseif ~internal.stats.isStrings(coefnames) % don't require cellstr
    error(message('stats:classreg:regr:modelutils:BadCoeffNames'));
elseif numel(coefnames) ~= ncoefs
    error(message('stats:classreg:regr:modelutils:BadCoeffNameLength'));
end

dfe = nobs - ncoefs;
tstat = coefs ./ se;
p = 2*tcdf(-abs(tstat),dfe);
t = table(coefs, se, tstat, p, ...
            'VariableNames',{'Estimate' 'SE' 'tStat' 'pValue'}, ...
            'RowNames',coefnames);
