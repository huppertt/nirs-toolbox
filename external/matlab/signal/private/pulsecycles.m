function [midRef,initCross,finalCross,nextCross] ...
   = pulsecycles(x, plotFlag, pulseTag, varargin)
%PULSECYCLES extract initial, final, and next pulse transitions
%
%   This function is for internal use only. It may be removed in the future.

%   Copyright 2011-2012 The MathWorks, Inc.

try
  % check if needs transposition
  needsTranspose = isrow(x);
  
  % extract leading numeric arguments
  [x, t, n] = chktransargs(0, x, varargin{:});
  
  % extract trailing optional arguments
  [tol, mprl, stateLevs, pol] = chktransopts(x, {'midpct','polarity'}, ...
                                           varargin{n:end});
catch ex
  throwAsCaller(ex);
end

% extract (absolute) reference levels from percent reference levels
[lwrBnd, uprBnd, midRef] = pctreflev(stateLevs,tol,100-tol, mprl);

% extract the mid-crossings and polarities with the specified state
% boundaries and reference levels
[midCross, polarity] = signal.internal.getMidCross(x, t, uprBnd, lwrBnd, midRef);

% find the first transition that matches the specified input polarity, pol.
iStart = find(polarity==pol, 1, 'first');
idx = iStart:numel(polarity);

if nargout>3
  % extract cycles (pulses inclusive of next transition)
  initCross = midCross(idx(1:2:end-2));
  finalCross = midCross(idx(2:2:end-1));
  nextCross = midCross(idx(3:2:end));
  allCross = [midCross(idx(1:2:end)); finalCross];    
else
  % extract pulses (exclusive of next transition)
  initCross = midCross(idx(1:2:end-1));
  finalCross = midCross(idx(2:2:end));
  allCross = [initCross; finalCross];
end

if plotFlag
  % plot all pulses of interest
  if strcmp(pulseTag,'PulseWidth')
    xdata = [initCross initCross finalCross finalCross]';
  elseif strcmp(pulseTag,'PulseSeparation')
    xdata = [finalCross finalCross nextCross nextCross]';
  elseif strcmp(pulseTag,'PulsePeriod')
    xdata = [initCross(1:2:end) initCross(1:2:end) nextCross(1:2:end) nextCross(1:2:end)]';
  else
    xdata = [];
  end
  ydata = repmat(stateLevs([1 2 2 1])',1,size(xdata,2));
  
  if ~isempty(xdata)
    [~, l] = signal.internal.plotInitialize(['signal:internal:plotInitialize:PlotName' pulseTag], x, t, ...
                  xdata, ydata, ['signal:internal:plotInitialize:PlotPatch' pulseTag]);
  else
    [~, l] = signal.internal.plotInitialize(['signal:internal:plotInitialize:PlotName' pulseTag], x, t);
  end
  l = signal.internal.plotItem(l,'signal:internal:plotItem:LegendMidCross', ...
                               allCross, midRef * ones(size(allCross)),'kx ', ...
                               'MarkerSize',10,'LineWidth',0.8);
  signal.internal.plotFinalize(l, [t(1) t(end)], stateLevs, tol, [NaN mprl NaN]);
end

if needsTranspose
  midRef     = midRef';
  initCross  = initCross';
  finalCross = finalCross';
  if nargout > 3
    nextCross  = nextCross';
  end
end