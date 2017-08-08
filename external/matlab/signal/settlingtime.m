function [s,sLev,sInst] = settlingtime(x, varargin)
%SETTLINGTIME Settling time metrics of bilevel waveform transitions
%   S = SETTLINGTIME(X,D) returns the duration of time each transition
%   takes to cross the mid-reference level to the point where the signal
%   enters and remains within a 2% tolerance of the final state over the
%   duration, D, specified as a real scalar. The settling times are
%   expressed as a vector, S, whose length corresponds to the number
%   of detected transitions in the input signal. If for any transition the
%   requested duration is not present or an intervening transition is
%   detected, the function marks the corresponding element in S as NaN. The
%   signal's sample values are specified by a vector, X, whose sample
%   instants correspond to the index of each sample value. To determine the
%   transitions, the function estimates the state levels of the input
%   waveform by a histogram method and identifies all regions which cross
%   between a 2% tolerance of each state level.
% 
%   S = SETTLINGTIME(X,Fs,D) specifies the sample rate, Fs, as a positive
%   scalar, where the first sample instant corresponds to a time of zero.
% 
%   S = SETTLINGTIME(X,T,D) specifies the sample instants, T, as a vector
%   with the same number of elements as X.
%
%   S = SETTLINGTIME(...,'Tolerance',TOL) specifies the tolerance that
%   the initial and final levels of each transition must be within their
%   respective state levels. The amount, TOL, is a scalar expressed
%   as the percentage of the difference between the upper and lower state
%   levels. The default value is 2.0 (percent).
% 
%   S = SETTLINGTIME(...,'MidPercentReferenceLevel',L) specifies the mid
%   reference level as a percentage of the waveform amplitude, where 0%
%   corresponds to the lower state level, and 100% corresponds to the upper
%   state level. L is a real scalar. The default value for L is 50
%   (percent).
% 
%   S = SETTLINGTIME(...,'StateLevels',SL) specifies the levels to use for
%   the lower and upper state levels, SL, as a two-element real row vector
%   whose first and second elements correspond to the lower and upper state
%   levels of the input waveform.
% 
%   [S,SLEV,SINST] = SETTLINGTIME(...) returns vectors, SLEV and
%   SINST, whose elements correspond to the level and sample instant of the
%   settling point of each transition.
%
%   SETTLINGTIME(...) plots the signal and darkens the regions of each
%   transition where settling time is computed.  It marks the location of the
%   settling time of each transition as well as the mid crossings and their
%   associated reference level. The state levels and their associated
%   lower and upper boundaries (adjustable by the TOL parameter) are also
%   plotted.
%
%   % Example 1:
%   %   Plot settling time information of a 2.3V step waveform 
%   %   over a duration of 2ns.
%   load('transitionex.mat', 'x', 't');
%   settlingtime(x, t, 2e-6)
%
%   % Example 2:
%   %   Compute settling time information of a 2.3V step waveform 
%   %   sampled at 4 MHz over a duration of 2ns.
%   load('transitionex.mat', 'x', 't');
%   y = settlingtime(x, 4e6, 2e-6)
%
% See also:  STATELEVELS, RISETIME, OVERSHOOT, MIDCROSS.

%   Copyright 2011-2013 The MathWorks, Inc.

% check to see if output needs transposition
needsTranspose = isrow(x);

% extract leading numeric arguments with delay parameter
[x, t, n, tSeek] = chktransargs(1, x, varargin{:});

% extract trailing optional arguments
[tol, mprl, stateLevs] = chktransopts(x, {'midpct'}, varargin{n:end});

% compute state boundaries and mid-reference level
[lwrBnd, uprBnd, midRef] = pctreflev(stateLevs,tol,100-tol, mprl);

% compute the mid crossings
[tm.MiddleCross, tm.Polarity, iPre, iPost] = signal.internal.getMidCross(x, t, uprBnd, lwrBnd, midRef);

% compute settling metrics
sm = signal.internal.getSettling(x, t, stateLevs(2), stateLevs(1), tol, tSeek, ...
                                 tm, iPre, iPost, false);

s = sm.Duration;
sLev = sm.Level;
sInst = sm.Instant;

if nargout==0
  % plot results if no output specified
  xdata = [tm.MiddleCross tm.MiddleCross sm.Instant sm.Instant]';
  ydata = [midRef * ones(size(tm.MiddleCross)) sm.Level sm.Level midRef * ones(size(tm.MiddleCross))]';
  [~, l] = signal.internal.plotInitialize( ...
      'signal:internal:plotInitialize:PlotNameSettlingTime', x, t, ...
      xdata, ydata, 'signal:internal:plotInitialize:PlotPatchSettlingTime');
  l = signal.internal.plotItem(l,'signal:internal:plotItem:LegendMidCross', ...
                               tm.MiddleCross, midRef * ones(size(s)),'kx ', ...
                               'MarkerSize',10,'LineWidth',0.8);
  l = signal.internal.plotItem(l,'signal:internal:plotItem:LegendSettlingPoint', ...
                               sInst, sLev,'ko ');
  signal.internal.plotFinalize(l, [t(1) t(end)], stateLevs, tol, [NaN mprl NaN]);
end

if needsTranspose
  s = s';
  sLev = sLev';
  sInst = sInst';
end