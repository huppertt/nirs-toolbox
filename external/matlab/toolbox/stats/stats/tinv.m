function x = tinv(p,v)
%TINV   Inverse of Student's T cumulative distribution function (cdf).
%   X=TINV(P,V) returns the inverse of Student's T cdf with V degrees 
%   of freedom, at the values in P.
%
%   The size of X is the common size of P and V. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   See also TCDF, TPDF, TRND, TSTAT, ICDF.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.6.2

%   Copyright 1993-2014 The MathWorks, Inc.


if nargin < 2, 
    error(message('stats:tinv:TooFewInputs')); 
end

[errorcode,p,v] = distchck(2,p,v);

if errorcode > 0
    error(message('stats:tinv:InputSizeMismatch'));
end

% Initialize Y to zero, or NaN for invalid d.f.
x = NaN(size(p),'like',internal.stats.dominantType(p,v)); % Will be single if either p or v is

% The inverse cdf of 0 is -Inf, and the inverse cdf of 1 is Inf.
x(p==0 & v > 0) = -Inf;
x(p==1 & v > 0) = Inf;

k0 = (0<p & p<1) & (v > 0);

% Invert the Cauchy distribution explicitly
k = find(k0 & (v == 1));
if any(k)
  x(k) = tan(pi * (p(k) - 0.5));
end

% For small d.f., call betaincinv which uses Newton's method
k = find(k0 & (v < 1000) & (v~=1));
if any(k)
    q = p(k) - .5;
    df = v(k);
    t = (abs(q) < .25);
    z = zeros(size(q), 'like', x);
    oneminusz = zeros(size(q), 'like', x);
    if any(t)
        % for z close to 1, compute 1-z directly to avoid roundoff
        oneminusz(t) = betaincinv(2.*abs(q(t)),0.5,df(t)/2,'lower');
        z(t) = 1 - oneminusz(t);
    end
    t = ~t; % (abs(q) >= .25);
    if any(t)
        z(t) = betaincinv(2.*abs(q(t)),df(t)/2,0.5,'upper');
        oneminusz(t) = 1 - z(t);
    end
    x(k) = sign(q) .* sqrt(df .* (oneminusz./z));
end

% For large d.f., use Abramowitz & Stegun formula 26.7.5
% k = find(p>0 & p<1 & ~isnan(x) & v >= 1000);
k = find(k0 & (v >= 1000));
if any(k)
   xn = norminv(p(k));
   df = v(k);
   x(k) = xn + (xn.^3+xn)./(4*df) + ...
           (5*xn.^5+16.*xn.^3+3*xn)./(96*df.^2) + ...
           (3*xn.^7+19*xn.^5+17*xn.^3-15*xn)./(384*df.^3) +...
           (79*xn.^9+776*xn.^7+1482*xn.^5-1920*xn.^3-945*xn)./(92160*df.^4);
end
