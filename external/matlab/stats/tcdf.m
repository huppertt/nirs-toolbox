function p = tcdf(x,v,uflag)
%TCDF   Student's T cumulative distribution function (cdf).
%   P = TCDF(X,V) computes the cdf for Student's T distribution
%   with V degrees of freedom, at the values in X.
%
%   The size of P is the common size of X and V. A scalar input
%   functions as a constant matrix of the same size as the other input.
%
%   P = TCDF(X,V,'upper') computes the upper tail probability of the
%   Student's T distribution with V degrees of freedom, at the values in X.
%
%   See also TINV, TPDF, TRND, TSTAT, CDF.

%   References:
%      [1] M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.7.
%      [2] L. Devroye, "Non-Uniform Random Variate Generation",
%      Springer-Verlag, 1986
%      [3] E. Kreyszig, "Introductory Mathematical Statistics",
%      John Wiley, 1970, Section 10.3, pages 144-146.

%   Copyright 1993-2013 The MathWorks, Inc.


normcutoff = 1e7;
if nargin < 2,
    error(message('stats:tcdf:TooFewInputs'));
end

[errorcode x v] = distchck(2,x,v);

if errorcode > 0
    error(message('stats:tcdf:InputSizeMismatch'));
end

if nargin>2
    % Compute upper tail
    if ~strcmpi(uflag,'upper')
        error(message('stats:cdf:UpperTailProblem'));
    end
    x = -x;
    uflag = [];
end

% Initialize P.
if isa(x,'single') || isa(v,'single')
    p = NaN(size(x),'single');
else
    p = NaN(size(x));
end

nans = (isnan(x) | ~(0<v)); %  v==NaN ==> (0<v)==false
cauchy = (v == 1);
normal = (v > normcutoff);

% General case: first compute F(-|x|) < .5, the lower tail.
%
% See Abramowitz and Stegun, formulas and 26.7.1/26.5.27 and 26.5.2
general = ~(cauchy | normal | nans);
xsq = x.^2;
% For small v, form v/(v+x^2) to maintain precision
t = (v < xsq) & general;
if any(t(:))
    p(t) = betainc(v(t) ./ (v(t) + xsq(t)), v(t)/2, 0.5, 'lower') / 2;
end

% For large v, form x^2/(v+x^2) to maintain precision
t = (v >= xsq) & general;
if any(t(:))
    p(t) = betainc(xsq(t) ./ (v(t) + xsq(t)), 0.5, v(t)/2, 'upper') / 2;
end

% For x > 0, F(x) = 1 - F(-|x|).
xpos = (x > 0);
p(xpos) = 1 - p(xpos); % p < .5, cancellation not a problem

% Special case for Cauchy distribution.  See Devroye pages 29 and 450.
% Note that instead of the usual Cauchy formula (atan x)/pi + 0.5, 
% we use acot(-x), which is equivalent and avoids roundoff error.
p(cauchy)  = xpos(cauchy) + acot(-x(cauchy))/pi; 

% Normal Approximation for very large nu.
p(normal) = normcdf(x(normal));

% Make the result exact for the median.
p(x == 0 & ~nans) = 0.5;
