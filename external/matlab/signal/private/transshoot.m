function [s,sLev,sInst] = transshoot(x, dirFlag, plotFlag, varargin)
%TRANSSHOOT Extract overshoot/undershoots of bilevel waveform transitions
% If dirFlag is +1 extract overshoots
% If dirFlag is -1 extract undershoots
%
% If plotFlag is true, plot and annotate the specified aberration.
%
%   This function is for internal use only. It may be removed in the future.

%   Copyright 2011-2012 The MathWorks, Inc.

try
  % check to see if output needs transposition
  needsTranspose = isrow(x);
  
  % extract leading numeric arguments
  [x, t, n] = chktransargs(0, x, varargin{:});

  % extract trailing optional arguments
  [tol, prl, stateLevs, rgn, seekFactor] = ...
    chktransopts(x, {'pctreflevels', 'region'}, varargin{n:end});
catch ex
  throwAsCaller(ex)
end

% compute state boundaries and reference levels
% (set mid reference level to mean of upper and lower reference level)
[lwrBnd, uprBnd, lwrRef, midRef, uprRef] = ...
    pctreflev(stateLevs, tol, 100-tol, prl(1), mean(prl), prl(2));

% extract the transitions
[tm, iPre, iPost] = signal.internal.getTransitions(x, t, uprBnd, lwrBnd, ...
                                                   uprRef, midRef, lwrRef);
% extract the metrics from the pre- or post-transition aberration region
if strcmpi(rgn,'preshoot')
  pm = signal.internal.getPreshoots(x, t, stateLevs(2), stateLevs(1), ...
                 uprBnd, lwrBnd, seekFactor, tm, iPre, iPost);
  msgId = 'signal:internal:plotItem:LegendPre';
else
  pm = signal.internal.getPostshoots(x, t, stateLevs(2), stateLevs(1),  ...
                 uprBnd, lwrBnd, seekFactor, tm, iPre, iPost, 0);
  msgId = 'signal:internal:plotItem:LegendPost';
end

if dirFlag>0
  s = pm.Overshoot;
  sLev = pm.OvershootLevel;
  sInst = pm.OvershootInstant;
  msgIdTag = 'Overshoot';
  sMarker = 'kv ';
else
  s = pm.Undershoot;
  sLev = pm.UndershootLevel;
  sInst = pm.UndershootInstant;
  msgIdTag = 'Undershoot';
  sMarker = 'k^ ';
end

if plotFlag
  % plot results
  [~, l] = signal.internal.plotInitialize( ...
        ['signal:internal:plotInitialize:PlotName' msgIdTag], x, t);
  l = signal.internal.plotItem(l,'signal:internal:plotItem:LegendUpperCross', ...
                               tm.UpperCross, uprRef * ones(size(s)),'rx ', ...
                               'MarkerSize',10,'LineWidth',0.8);
  l = signal.internal.plotItem(l,'signal:internal:plotItem:LegendLowerCross', ...
                               tm.LowerCross, lwrRef * ones(size(s)),'gx ', ...
                               'MarkerSize',10,'LineWidth',0.8);
  l = signal.internal.plotItem(l,[msgId msgIdTag], sInst, sLev, sMarker);
  signal.internal.plotFinalize(l, [t(1) t(end)], stateLevs, tol, [prl(1) NaN prl(2)]);
end

if needsTranspose
  s     = s';
  sLev  = sLev';
  sInst = sInst';
end