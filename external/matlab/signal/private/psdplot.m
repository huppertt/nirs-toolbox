function [hLine, xscale] = psdplot(Pxx, F, rbw, esttype, status)
%PSDPLOT Helper function for plotting power and psd estimates.

%   Copyright 2013 The MathWorks, Inc.

hAxes = newplot;

if nargin>4 && status.normF
  F = F/pi;
  [F, xscale, xunits] = engunits(F);
  xscale = xscale/pi;
  normFreqStr = getString(message('signal:sigtools:getfrequnitstrs:NormalizedFrequency'));
  xlbl = [normFreqStr '  (\times\pi ' xunits 'rad/sample)'];
else
  [F, xscale, xunits] = engunits(F);  
  xlbl = getfreqlbl([xunits 'Hz']);
end


xlabel(hAxes, xlbl);


if strcmp(esttype,'power')
  H = 10*log10(rbw*Pxx);
  ylbl = getString(message('signal:dspdata:dspdata:PowerdB'));
else %'psd'
  H = 10*log10(Pxx);
  if nargin>4 && status.normF
    ylbl = getString(message('signal:dspdata:dspdata:PowerfrequencydBradsample'));
  else
    ylbl = getString(message('signal:dspdata:dspdata:PowerfrequencydBHz'));
  end
end

ylabel(hAxes, ylbl);

hLine = line(F, H, 'Parent', hAxes);
set(hAxes, 'XLim', [min(F) max(F)]);

% Ensure axes limits are properly cached for zoom/unzoom
resetplotview(hAxes,'SaveCurrentView');  

initdistgrid(hAxes);
uistack(hLine,'top');