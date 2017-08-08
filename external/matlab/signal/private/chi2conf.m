function c = chi2conf(conf,k)
%CHI2CONF Confidence interval using inverse of chi-square cdf
%   This is a helper function for Spectrum objects

%   C = CHI2CONF(CONF,K) calculates confidence interval C based on
%   confidence level CONF and K independent measurements.  C is a two
%   element vector.  
%
%   Reference:
%     Stephen Kay, "Modern Spectral Estimation, Theory & Application," 
%     Prentice Hall, 1988, pp 76, eqn 4.16. 

%   Copyright 2007-2013 The MathWorks, Inc.

narginchk(2,2);

% Ensure double precision arithmetic
conf = double(conf);
k = double(k);

v=2*k;
alfa = 1 - conf;
c=chi2inv([1-alfa/2 alfa/2],v);
c=v./c;

%--------------------------------------------------------------------------

function x = chi2inv(p,v)
%CHI2INV Inverse of the chi-square cumulative distribution function (cdf).
%   X = CHI2INV(P,V)  returns the inverse of the chi-square cdf with V  
%   degrees of freedom at the values in P. The chi-square cdf with V 
%   degrees of freedom, is the gamma cdf with parameters V/2 and 2.   
%
%   The size of X is the common size of P and V. A scalar input
%   functions as a constant matrix of the same size as the other input.   

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.4.
%      [2] E. Kreyszig, "Introductory Mathematical Statistics",
%      John Wiley, 1970, section 10.2 (page 144)


error(nargchk(2,2,nargin,'struct'));

[errorcode p v] = distchck(2,p,v);

if errorcode > 0
    error(message('signal:chi2conf:InvalidNonScalarDimensions'));
end

% Call the gamma inverse function. 
x = gaminv(p,v/2,2);

% Return NaN if the degrees of freedom is not a positive integer.
k = find(v < 0  |  round(v) ~= v);
if any(k)
    tmp  = NaN;
    x(k) = tmp(ones(size(k)));
end


function x = gaminv(p,a,b)
%GAMINV Inverse of the gamma cumulative distribution function (cdf).
%   X = GAMINV(P,A,B)  returns the inverse of the gamma cdf with  
%   parameters A and B, at the probabilities in P.
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   GAMINV uses Newton's method to converge to the solution.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 6.5.

%   B.A. Jones 1-12-93
%   Was: Revision: 1.2, Date: 1996/07/25 16:23:36

if nargin<3, 
    b=1;
end

[errorcode p a b] = distchck(3,p,a,b);

if errorcode > 0
    error(message('signal:chi2conf:InvalidDimensions'));
end

%   Initialize X to zero.
x = zeros(size(p));

k = find(p<0 | p>1 | a <= 0 | b <= 0);
if any(k),
    tmp = NaN;
    x(k) = tmp(ones(size(k)));
end

% The inverse cdf of 0 is 0, and the inverse cdf of 1 is 1.  
k0 = find(p == 0 & a > 0 & b > 0);
if any(k0),
    x(k0) = zeros(size(k0)); 
end

k1 = find(p == 1 & a > 0 & b > 0);
if any(k1), 
    tmp = Inf;
    x(k1) = tmp(ones(size(k1))); 
end

% Newton's Method
% Permit no more than count_limit iterations.
count_limit = 100;
count = 0;

k = find(p > 0  &  p < 1 & a > 0 & b > 0);
pk = p(k);

% Supply a starting guess for the iteration.
%   Use a method of moments fit to the lognormal distribution. 
mn = a(k) .* b(k);
v = mn .* b(k);
temp = log(v + mn .^ 2); 
mu = 2 * log(mn) - 0.5 * temp;
sigma = -2 * log(mn) + temp;
xk = exp(norminv(pk,mu,sigma));

h = ones(size(pk)); 

% Break out of the iteration loop for three reasons:
%  1) the last update is very small (compared to x)
%  2) the last update is very small (compared to sqrt(eps))
%  3) There are more than 100 iterations. This should NEVER happen. 

while(any(abs(h) > sqrt(eps)*abs(xk))  &  max(abs(h)) > sqrt(eps)    ...
                                 & count < count_limit), 
                                 
    count = count + 1;
    h = (gamcdf(xk,a(k),b(k)) - pk) ./ gampdf(xk,a(k),b(k));
    xnew = xk - h;
    % Make sure that the current guess stays greater than zero.
    % When Newton's Method suggests steps that lead to negative guesses
    % take a step 9/10ths of the way to zero:
    ksmall = find(xnew < 0);
    if any(ksmall),
        xnew(ksmall) = xk(ksmall) / 10;
        h = xk-xnew;
    end
    xk = xnew;
end


% Store the converged value in the correct place
x(k) = xk;

if count == count_limit, 
    fprintf(getString(message('signal:chi2conf:WarningGAMINVDidNotConverge')));
    str = getString(message('signal:chi2conf:TheLastStepWas'));
    outstr = sprintf([str,'  %13.8f'],h);
    fprintf(outstr);
end

function x = norminv(p,mu,sigma)
%NORMINV Inverse of the normal cumulative distribution function (cdf).
%   X = NORMINV(P,MU,SIGMA) finds the inverse of the normal cdf with
%   mean, MU, and standard deviation, SIGMA.
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   Default values for MU and SIGMA are 0 and 1 respectively.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 7.1.1 and 26.2.2

%   Was: Revision: 1.2, Date: 1996/07/25 16:23:36

if nargin < 3, 
    sigma = 1;
end

if nargin < 2;
    mu = 0;
end

[errorcode p mu sigma] = distchck(3,p,mu,sigma);

if errorcode > 0
    error(message('signal:chi2conf:InvalidNonScalarDimensions'));
end

% Allocate space for x.
x = zeros(size(p));

% Return NaN if the arguments are outside their respective limits.
k = find(sigma <= 0 | p < 0 | p > 1);
if any(k)
    tmp  = NaN;
    x(k) = tmp(ones(size(k))); 
end

% Put in the correct values when P is either 0 or 1.
k = find(p == 0);
if any(k)
    tmp  = Inf;
    x(k) = -tmp(ones(size(k)));
end

k = find(p == 1);
if any(k)
    tmp  = Inf;
    x(k) = tmp(ones(size(k))); 
end

% Compute the inverse function for the intermediate values.
k = find(p > 0  &  p < 1 & sigma > 0);
if any(k),
    x(k) = sqrt(2) * sigma(k) .* erfinv(2 * p(k) - 1) + mu(k);
end


function p = gamcdf(x,a,b)
%GAMCDF Gamma cumulative distribution function.
%   P = GAMCDF(X,A,B) returns the gamma cumulative distribution
%   function with parameters A and B at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   Some references refer to the gamma distribution with a single
%   parameter. This corresponds to the default of B = 1. 
%
%   GAMMAINC does computational work.

%   References:
%      [1]  L. Devroye, "Non-Uniform Random Variate Generation", 
%      Springer-Verlag, 1986. p. 401.
%      [2]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.32.

%   Was: Revision: 1.2, Date: 1996/07/25 16:23:36

error(nargchk(2,3,nargin,'struct'));

if nargin < 3, 
    b = 1; 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error(message('signal:chi2conf:InvalidNonScalarDimensions'));
end

%   Return NaN if the arguments are outside their respective limits.
k1 = find(a <= 0 | b <= 0);     
if any(k1)
    tmp   = NaN;
    p(k1) = tmp(ones(size(k1)));
end

% Initialize P to zero.
p = zeros(size(x));

k = find(x > 0 & ~(a <= 0 | b <= 0));
if any(k), 
    p(k) = gammainc(x(k) ./ b(k),a(k));
end

% Make sure that round-off errors never make P greater than 1.
k = find(p > 1);
if any(k)
    p(k) = ones(size(k));
end


function y = gampdf(x,a,b)
%GAMPDF Gamma probability density function.
%   Y = GAMPDF(X,A,B) returns the gamma probability density function 
%   with parameters A and B, at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   Some references refer to the gamma distribution with a single
%   parameter. This corresponds to the default of B = 1.

%   References:
%      [1]  L. Devroye, "Non-Uniform Random Variate Generation", 
%      Springer-Verlag, 1986, pages 401-402.

%   Was: Revision: 1.2, Date: 1996/07/25 16:23:36

error(nargchk(2,3,nargin,'struct'));

if nargin < 3, 
    b = 1; 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error(message('signal:chi2conf:InvalidNonScalarDimensions'));
end

% Initialize Y to zero.
y = zeros(size(x));

%   Return NaN if the arguments are outside their respective limits.
k1 = find(a <= 0 | b <= 0);     
if any(k1)
    tmp = NaN;
    y(k1) = tmp(ones(size(k1)));
end

k=find(x > 0 & ~(a <= 0 | b <= 0));
if any(k)
    y(k) = (a(k) - 1) .* log(x(k)) - (x(k) ./ b(k)) - gammaln(a(k)) - a(k) .* log(b(k));
    y(k) = exp(y(k));
end
k1 = find(x == 0 & a < 1);
if any(k1)
  tmp = Inf;
  y(k1) = tmp(ones(size(k1)));
end
k2 = find(x == 0 & a == 1);
if any(k2)
  y(k2) = (1./b(k2));
end

function [errorcode,out1,out2,out3,out4] = distchck(nparms,arg1,arg2,arg3,arg4)
%DISTCHCK Checks the argument list for the probability functions.

%   B.A. Jones  1-22-93
%   Was: Revision: 1.2, Date: 1996/07/25 16:23:36

errorcode = 0;

if nparms == 1
    out1 = arg1;
    return;
end
    
if nparms == 2
    [r1 c1] = size(arg1);
    [r2 c2] = size(arg2);
    scalararg1 = (prod(size(arg1)) == 1);
    scalararg2 = (prod(size(arg2)) == 1);
    if ~scalararg1 & ~scalararg2
        if r1 ~= r2 | c1 ~= c2
            errorcode = 1;
            return;         
        end     
    end
    if scalararg1
        out1 = arg1(ones(r2,1),ones(c2,1));
    else
        out1 = arg1;
    end
    if scalararg2
        out2 = arg2(ones(r1,1),ones(c1,1));
    else
        out2 = arg2;
    end
end
    
if nparms == 3
    [r1 c1] = size(arg1);
    [r2 c2] = size(arg2);
    [r3 c3] = size(arg3);
    scalararg1 = (prod(size(arg1)) == 1);
    scalararg2 = (prod(size(arg2)) == 1);
    scalararg3 = (prod(size(arg3)) == 1);

    if ~scalararg1 & ~scalararg2
        if r1 ~= r2 | c1 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg1 & ~scalararg3
        if r1 ~= r3 | c1 ~= c3
            errorcode = 1;
            return;                 
        end
    end

    if ~scalararg3 & ~scalararg2
        if r3 ~= r2 | c3 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg1
      out1 = arg1;
    end
    if ~scalararg2
      out2 = arg2;
    end
    if ~scalararg3
      out3 = arg3;
    end
    rows = max([r1 r2 r3]);
   columns = max([c1 c2 c3]);  
       
    if scalararg1
        out1 = arg1(ones(rows,1),ones(columns,1));
    end
   if scalararg2
        out2 = arg2(ones(rows,1),ones(columns,1));
   end
   if scalararg3
       out3 = arg3(ones(rows,1),ones(columns,1));
   end
     out4 =[];
    
end

if nparms == 4
    [r1 c1] = size(arg1);
    [r2 c2] = size(arg2);
    [r3 c3] = size(arg3);
    [r4 c4] = size(arg4);
    scalararg1 = (prod(size(arg1)) == 1);
    scalararg2 = (prod(size(arg2)) == 1);
    scalararg3 = (prod(size(arg3)) == 1);
    scalararg4 = (prod(size(arg4)) == 1);

    if ~scalararg1 & ~scalararg2
        if r1 ~= r2 | c1 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg1 & ~scalararg3
        if r1 ~= r3 | c1 ~= c3
            errorcode = 1;
            return;                 
        end
    end

    if ~scalararg1 & ~scalararg4
        if r1 ~= r4 | c1 ~= c4
            errorcode = 1;
            return;                 
        end
    end

    if ~scalararg3 & ~scalararg2
        if r3 ~= r2 | c3 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg4 & ~scalararg2
        if r4 ~= r2 | c4 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg3 & ~scalararg4
        if r3 ~= r4 | c3 ~= c4
            errorcode = 1;
            return;         
        end
    end


    if ~scalararg1
       out1 = arg1;
    end
    if ~scalararg2
       out2 = arg2;
    end
    if ~scalararg3
      out3 = arg3;
    end
    if ~scalararg4
      out4 = arg4;
    end
 
   rows = max([r1 r2 r3 r4]);
   columns = max([c1 c2 c3 c4]);     
    if scalararg1
       out1 = arg1(ones(rows,1),ones(columns,1));
   end
   if scalararg2
       out2 = arg2(ones(rows,1),ones(columns,1));
   end
   if scalararg3
       out3 = arg3(ones(rows,1),ones(columns,1));
   end
   if scalararg4
       out4 = arg4(ones(rows,1),ones(columns,1));
   end
end



% [EOF]
