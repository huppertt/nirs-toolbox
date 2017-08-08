function [y,midRef] = midcross(x, varargin)
%MIDCROSS Mid reference level crossing of bilevel waveform transitions
%   C = MIDCROSS(X) returns linearly interpolated time instants where each
%   transition of the input signal crosses the 50% reference level. The
%   signal's sample values are specified by a vector, X, whose sample
%   instants correspond to the index of each sample value. To determine the
%   transitions, the function estimates the state levels of the input
%   waveform by a histogram method and identifies all regions which cross
%   between a 2% tolerance of each state level.
%
%   C = MIDCROSS(X,Fs) specifies the sample rate, Fs, as a positive
%   scalar, where the first sample instant corresponds to a time of zero.
%
%   C = MIDCROSS(X,T) specifies the sample instants, T, as a vector with
%   the same number of elements as X.
%
%   C = MIDCROSS(...,'Tolerance',TOL) specifies the tolerance that the
%   initial and final levels of each transition must be within their
%   respective state levels. The amount, TOL, is a real scalar expressed as
%   the percentage of the difference between the upper and lower state
%   levels. The default value is 2.0 (percent).
%
%   C = MIDCROSS(...,'MidPercentReferenceLevel',L) specifies the mid
%   reference level as a percentage of the waveform amplitude, where 0%
%   corresponds to the lower state level, and 100% corresponds to the upper
%   state level. L is a real scalar. The default value for L is 50
%   (percent).
%
%   C = MIDCROSS(...,'StateLevels',SL) specifies the levels to use for the
%   lower and upper state levels, SL, as a two-element real row vector
%   whose first and second elements correspond to the lower and upper state
%   levels of the input waveform.
%
%   [C, MIDLEV] = MIDCROSS(...) returns the scalar, MIDLEV, that
%   corresponds to the mid-reference level.
%
%   MIDCROSS(...) plots the signal and marks the location of the mid
%   crossings, and their associated reference level. The state levels and
%   their associated lower and upper boundaries (adjustable by the TOL
%   parameter) are also plotted.
%
%   % Example 1:
%   %   Plot the mid crossing information of a 2.3V step waveform 
%   load('transitionex.mat', 'x', 't');
%   midcross(x, t)
%
%   % Example 2:
%   %   Compute the mid crossing information of a 2.3V step waveform 
%   %   sampled at 4 MHz
%   load('transitionex.mat', 'x', 't');
%   y = midcross(x, 4e6)
%
% See also:  STATELEVELS, RISETIME, OVERSHOOT, SETTLINGTIME, PULSEWIDTH.

%   Copyright 2011-2013 The MathWorks, Inc.

% check to see if output needs transposition
needsTranspose = isrow(x);

% extract leading numeric arguments
[x, t, n] = chktransargs(0, x, varargin{:});

% extract trailing optional arguments
[tol, mprl, stateLevs] = chktransopts(x, {'midpct'}, varargin{n:end});

% compute state boundaries and mid reference level
[lwrBnd, uprBnd, midRef] = pctreflev(stateLevs, tol, 100-tol, mprl);
    
% extract the mid-crossings
y = signal.internal.getMidCross(x, t, uprBnd, lwrBnd, midRef);

if nargout==0
  % plot results if no output specified
  [~, l] = signal.internal.plotInitialize('signal:internal:plotInitialize:PlotNameMidCross', x, t);
  l = signal.internal.plotItem(l,'signal:internal:plotItem:LegendMidCross', ...
                               y, midRef * ones(size(y)),'kx ', ...
                               'MarkerSize',10,'LineWidth',0.8);
  signal.internal.plotFinalize(l, [t(1) t(end)], stateLevs, tol, [NaN mprl NaN]);
end

if needsTranspose
  y = y';
end
