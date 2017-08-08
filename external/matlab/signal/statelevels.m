function [levels,histogram,bins] = statelevels(x, nbins, method, bounds)
%STATELEVELS State level estimation of bilevel waveforms via histogram
%   LEVELS = STATELEVELS(X) returns the lower and upper state levels in the
%   real vector, X.  The levels are estimated by computing a histogram
%   which contains 100 bins and whose lower and upper bounds correspond to
%   the minimum and maximum value found in X.  The computed histogram is
%   divided into two equal sized regions between the first and last
%   non-zero bin. The mode of each region of the histogram is returned as
%   the first and second element of a two-element row vector, LEVELS.
%   
%   LEVELS = STATELEVELS(X,NBINS) specifies the number of histogram bins to 
%   use in the histogram as a positive scalar.
%   
%   LEVELS = STATELEVELS(X,NBINS,METHOD) performs state level estimation by
%   computing the specified METHOD as one of 'mean' | 'mode', where the
%   default value is 'mode'.  If 'mean' or 'mode' is specified, the mean or
%   mode of each of the two regions of the computed histogram is returned
%   as the first and second element of LEVELS, respectively.
%
%   LEVELS = STATELEVELS(X,NBINS,METHOD,BOUNDS) specifies the lower and
%   upper bounds of the histogram as the first and second element of the
%   real, two-element row vector BOUNDS.  Extreme values of X which lie
%   outside these bounds are ignored when computing the histogram.
%
%   [LEVELS,HISTOGRAM] = STATELEVELS(...) returns a vector containing the
%   computed histogram, HISTOGRAM.
%
%   [LEVELS,HISTOGRAM,BINLEVELS] = STATELEVELS(...) returns a vector
%   of levels, BINLEVELS, that correspond to the center of each bin in the
%   histogram.
%
%   STATELEVELS(...) displays a plot of the signal and the corresponding
%   computed histogram.
%   
%   % Example 1:
%   %   Plot the state levels of a 2.3V under-damped noisy clock.
%   load('clockex.mat', 'x', 't');
%   statelevels(x)
%
%   % Example 2:
%   %   Compute the state levels of a 2.3V under-damped noisy clock.
%   load('clockex.mat', 'x', 't');
%   levels = statelevels(x)
%
%   See also MIDCROSS RISETIME OVERSHOOT PULSEWIDTH.

%   Copyright 2011-2012 The MathWorks, Inc.

validateattributes(x,{'double'},{'real','finite','vector'}, ...
                   'statelevels','X');
if numel(x) < 2
  error(message('signal:statelevels:MustBeMultiElementVector','X'));
end

if nargin > 1
  validateattributes(nbins,{'double'},{'real','finite','scalar','>',1}, ...
                     'statelevels','NBINS');
else
  nbins = 100;
end

if nargin > 2
  method = validatestring(method,{'mean','mode'},'statelevels','METHOD');
  method = ['Histogram ' method];
else
  method = 'Histogram mode';
end

if nargin > 3
  validateattributes(bounds,{'double'},{'real','finite','size',[1 2]}, ...
                     'statelevels','BOUNDS');
  lower = bounds(1);
  upper = bounds(2);
  if lower >= upper
    error(message('signal:statelevels:MustBeStrictlyIncreasing','BOUNDS'));
  end
else
  lower = min(x);
  upper = max(x);
end

histogram = signal.internal.getHistogram(x(:), nbins, lower, upper);
levels = signal.internal.getLevelsByHistogram(histogram, lower, upper, method);
bins = lower + ((1:nbins) - 0.5)' * (upper - lower) / nbins;

if nargout == 0
  signal.internal.plotStateLevels(x, levels, false, method, lower, upper, histogram);
end    

if isrow(x)
  histogram = histogram.';
  bins = bins.';
end