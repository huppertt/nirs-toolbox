function [hAxes, F, fscale, colorOrder] = initdistplot(plotType, F)
%INITDISTPLOT initialize a plot for distortion annotation
%
%   This function is for internal use only. It may be removed in the future.
   
%   Copyright 2013 The MathWorks, Inc.

hAxes = newplot;

% set up x axis
[F, fscale, xunits] = engunits(F);
xlbl = getfreqlbl([xunits 'Hz']);
xlabel(hAxes, xlbl);

% bound plot by x limit
set(hAxes, 'XLim', [min(F) max(F)]);

% set up y axis
if strcmp(plotType,'power')
  ylbl = getString(message('signal:dspdata:dspdata:PowerdB'));
else %'psd'
  ylbl = getString(message('signal:dspdata:dspdata:PowerfrequencydBHz'));
end
ylabel(hAxes, ylbl);

% use the default color order (but set third color to black instead of red)
colorOrder = get(0,'DefaultAxesColorOrder');
colorOrder(3,:) = [0 0 0];