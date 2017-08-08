function CI = arconfinterval(CL, p, Pxx, x)
%ARCONFINTERVAL  Confidence Interval for AR spectrum estimation methods.
%   CI = ARCONFINTERVAL(CL,ORDER,PXX,X) calculates the confidence
%   interval CI for spectrum estimate PXX based on confidence level CL. 
%   X is the data used for computing the spectrum estimate PXX.
% 
%   Reference : Steven M.Kay, "Modern spectral Estimation",
%   Prentice Hall, 1988, Chapter 6, pp 194-195
% 
%   Copyright 2012-2014 The MathWorks, Inc.

% Compute intervals using double precision arithmetic
CL = double(CL);
p = double(p);
Pxx = double(Pxx);
x = double(x);

alfa = 1-CL;
normval = norminverse(1-alfa/2,0,1);
if isvector(x)
  N = length(x);
else
  N = size(x,1);
end

if( N/(2*p) > normval^2)
    beta = sqrt(2*p/N)* normval;
    % do upper bound first to completely pre-allocate CI.
    CI(:,2:2:2*size(Pxx,2)) = Pxx*(1+beta);
    CI(:,1:2:2*size(Pxx,2)) = Pxx*(1-beta);
else
    warning(message('signal:spectrum:abstractar:confinterval:InsufficientData'));
    CI = [];
end


%--------------------------------------------------------------------------
function [x] = norminverse(p,mu,sigma)
%NORMINVERSE Inverse of the normal cumulative distribution function (cdf).
%   X = NORMINVERSE(P,MU,SIGMA) returns the inverse cdf for the normal
%   distribution with mean MU and standard deviation SIGMA, evaluated at
%   the value in P.  
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%
%   References:
%      [1] Abramowitz, M. and Stegun, I.A. (1964) Handbook of Mathematical
%          Functions, Dover, New York, 1046pp., sections 7.1, 26.2.
%      [2] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley, 170pp.

if(p < 0 || 1 < p)
     error(message('signal:spectrum:abstractar:confinterval:InvalidValue', 'P'));
end

x0 = -sqrt(2).*erfcinv(2*p);
x = sigma*x0 + mu;
