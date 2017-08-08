function [z,mu,sigma] = wnanzscore(x,w,biased)
%WNANZSCORE Weighted z-score ignoring NaN's.
%   Z=WNANZSCORE(X,W) returns a centered, scaled version of X, the same
%   size as X. For vector input X, Z is the vector of z-scores
%   (X-MU)./SIGMA, where MU is the weighted mean and SIGMA is the weighted
%   standard deviation computed ignoring NaN's. For matrix X, z-scores are
%   computed using the weighted mean and standard deviation along each
%   column of X. Pass W as a vector with one element per row of X.
%
%   [Z,MU,SIGMA]=WNANZSCORE(X,W) also returns the weighted mean and
%   standard deviation.
%
%   Z=WNANZSCORE(X,W,1) normalizes X using biased (maximum likelihood)
%   estimates of MU and SIGMA. Z=WNANZSCORE(X,W,0) is the same as
%   Z=WNANZSCORE(X,W). Note that WNANVAR uses the opposite convention.
%
%   See also classreg.learning.internal.wnanmean,
%   classreg.learning.internal.wnanvar.

% [] is a special case for std and mean, just handle it out here.
if isequal(x,[]), z = []; return; end

% Default is unbiased std
if nargin<3
    biased = 0;
end

% Compute X's mean and sd, and standardize it
mu = classreg.learning.internal.wnanmean(x,w);
sigma = sqrt(classreg.learning.internal.wnanvar(x,w,~biased));    % wnanvar uses the opposite semantics for its 'flag'.
sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = bsxfun(@minus,x, mu);
z = bsxfun(@rdivide, z, sigma0);
end