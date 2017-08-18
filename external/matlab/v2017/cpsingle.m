function [icp, residue] = cpsingle(y, statistic, Lmin)
%CPSINGLE single best changepoint between two equal distributions
%   This file is for internal use only and may be removed in a future
%   release.
%
%   References:
%   * Tony F. Chan, Gene H. Golub, Randall J. LeVeque, Algorithms for
%     computing the sample variance: analysis and recommendations"  The
%     American Statistician.  Vol 37, No. 3 (Aug., 1983) pp. 242-247.
%   * Philippe Pierre Pebay. "Formulas for robust one-pass parallel
%     computation of covariances and arbitrary-order statistical moments."
%     SAND2008-6212.  2008. Published through SciTech Connect.

%   Copyright 2015 The MathWorks, Inc.


% farm out to correct algorithm
if strcmp(statistic,'mean')
  [fwd, rev] = meanResidual(y);
elseif strcmp(statistic,'rms')
  [fwd, rev] = rmsResidual(y);
elseif strcmp(statistic,'std')
  [fwd, rev] = stdResidual(y);
elseif strcmp(statistic,'linear')
  [fwd, rev] = linearResidual(y);
end

[residue, icp] = min(fwd(Lmin:end-Lmin)+rev(Lmin+1:end-Lmin+1));
icp = icp + Lmin;

% ------------------------------------------------------------------------
function [fwd,rev] = meanResidual(y)
fwd = zeros(size(y));
rev = zeros(size(y));
n = numel(y);

ymean = 0;
Syy = 0;
for ix=1:n
  ydelta = y(ix) - ymean;
  npoints = ix;
  ymean = ymean + ydelta./npoints;
  Syy = Syy + ydelta * (y(ix) - ymean);
  fwd(ix) = Syy;
end

ymean = 0;
Syy = 0;
for ix=n:-1:1
  ydelta = y(ix) - ymean;
  npoints = n-ix+1;
  ymean = ymean + ydelta./npoints;
  Syy = Syy + ydelta * (y(ix) - ymean);
  rev(ix) = Syy;
end

% ------------------------------------------------------------------------
function [fwd,rev] = rmsResidual(y)
fwd = zeros(size(y));
rev = zeros(size(y));
n = numel(y);

logrealmin = log(realmin);

Syy = 0;
for ix=1:n
  npoints = ix;
  Syy = Syy + y(ix).^2;
  fwd(ix) = npoints*max(logrealmin,log(Syy./npoints));
end

Syy = 0;
for ix=n:-1:1
  npoints = n+1-ix;
  Syy = Syy + y(ix).^2;
  rev(ix) = npoints*max(logrealmin,log(Syy./npoints));
end

% ------------------------------------------------------------------------
function [fwd,rev] = stdResidual(y)
fwd = zeros(size(y));
rev = zeros(size(y));
n = numel(y);

logrealmin = log(realmin);

ymean = 0;
Syy = 0;
for ix=1:n
  npoints = ix;
  ydelta = y(ix) - ymean;
  ymean = ymean + ydelta./npoints;
  Syy = Syy + ydelta * (y(ix) - ymean);
  fwd(ix) = npoints*max(logrealmin,log(Syy./npoints));
end

ymean = 0;
Syy = 0;
for ix=n:-1:1
  ydelta = y(ix) - ymean;
  npoints = n+1-ix;
  ymean = ymean + ydelta./npoints;
  Syy = Syy + ydelta * (y(ix) - ymean);
  rev(ix) = npoints*max(logrealmin,log(Syy./npoints));
end

% ------------------------------------------------------------------------
function [fwd,rev] = linearResidual(y)
fwd = zeros(size(y));
rev = zeros(size(y));

xmean = 0;
ymean = 0;
Sxx = 0;
Syy = 0;
Sxy = 0;
SxxSSE = 0;

for ix=1:numel(y)
  npoints = ix;
  ydelta = y(ix) - ymean;
  xdelta = ix - xmean;
  ymean = ymean + ydelta/npoints;
  xmean = xmean + xdelta/npoints;
  dSyy = ydelta .* (y(ix) - ymean);
  dSxx = xdelta .* (ix - xmean);
  dSxy = xdelta .* ydelta .* (npoints - 1) ./ npoints;
  Syy = Syy + dSyy;
  dSxxSSE = dSxx.*Syy + dSyy.*Sxx - dSxy.*(2*Sxy+dSxy);
  Sxx = Sxx + dSxx;
  Sxy = Sxy + dSxy;
  SxxSSE = SxxSSE + dSxxSSE;
  fwd(ix) = SxxSSE ./ Sxx;
end

xmean = 0;
ymean = 0;
Sxx = 0;
Syy = 0;
Sxy = 0;
n = numel(y);
for ix=n:-1:1
  npoints = n+1-ix;
  ydelta = y(ix) - ymean;
  xdelta = ix - xmean;
  ymean = ymean + ydelta/npoints;
  xmean = xmean + xdelta/npoints;
  Syy = Syy + ydelta * (y(ix) - ymean);
  Sxx = Sxx + xdelta * (ix - xmean);
  Sxy = Sxy + xdelta .* ydelta .* (npoints - 1) ./ npoints;
  rev(ix) = Syy - Sxy.^2 ./ Sxx;
end
