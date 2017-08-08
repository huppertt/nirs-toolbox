function [d,lwrCross,uprCross,lwrRef,uprRef] ...
    = transdurs(x, polarityFlag, plotFlag, varargin)
%TRANSDURS Extract durations of bilevel waveform transitions
%   If polarityFlag is +1 return duration metrics for positive transitions
%   If polarityFlag is -1 return duration metrics for negative transitions
%   If polarityFlag is  0 return duration metrics for all transitions
%
%   If plotFlag is true, plot the locations of the upper and lower
%   crossings.
%
%   This function is for internal use only. It may be removed in the future.

%   Copyright 2011-2012 The MathWorks, Inc.

try 
  % check to see if input is a row vector
  needsTranspose = isrow(x);
  
  % extract leading numeric arguments
  [x, t, n] = chktransargs(0, x, varargin{:});

  % extract trailing optional arguments
  [tol, prl, stateLevs] = chktransopts(x, {'PctRefLevels'}, varargin{n:end});
catch ex
  throwAsCaller(ex);
end
  
% compute state boundaries and reference levels
% (set mid reference level to mean of upper and lower reference level)
[lwrBnd, uprBnd, lwrRef, midRef, uprRef] = ...
  pctreflev(stateLevs, tol, 100-tol, prl(1), mean(prl), prl(2));

% extract the transitions
tm = signal.internal.getTransitions(x, t, uprBnd, lwrBnd, uprRef, midRef, lwrRef);

% find the desired polarities
if polarityFlag<0
  idx = find(tm.Polarity<0);
  msgIdTag = 'FallTime';
elseif polarityFlag>0
  idx = find(tm.Polarity>0);
  msgIdTag = 'RiseTime';
else
  idx = 1:numel(tm.Polarity);
  msgIdTag = 'SlewRate';
end

d = tm.Duration(idx);
lwrCross = tm.LowerCross(idx);
uprCross = tm.UpperCross(idx);

if plotFlag
  xdata = [lwrCross lwrCross uprCross uprCross]';
  ydata = repmat([lwrRef uprRef uprRef lwrRef]',1,numel(lwrCross));  
  [~, l] = signal.internal.plotInitialize( ...
            ['signal:internal:plotInitialize:PlotName' msgIdTag], x, t, ...
            xdata, ydata, ['signal:internal:plotInitialize:PlotPatch' msgIdTag]);
  l = signal.internal.plotItem(l,'signal:internal:plotItem:LegendUpperCross', ...
                               uprCross, uprRef * ones(size(idx)),'rx ', ...
                               'MarkerSize',10,'LineWidth',0.8);
  l = signal.internal.plotItem(l,'signal:internal:plotItem:LegendLowerCross', ...
                               lwrCross, lwrRef * ones(size(idx)),'gx ', ...
                               'MarkerSize',10,'LineWidth',0.8);
  signal.internal.plotFinalize(l, [t(1) t(end)], stateLevs, tol, [prl(1) NaN prl(2)]);
end

if needsTranspose
  d = d';
  lwrCross = lwrCross';
  uprCross = uprCross';
end
  