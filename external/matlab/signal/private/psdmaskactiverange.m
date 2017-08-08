function psdmaskactiverange(hAxes, xscale, xLim, yLim, xRange)
%PSDMASKACTIVERANGE Helper function for ranged power and psd estimates
%   Masks the active range via shading with transparent patches

%   Copyright 2014 The MathWorks, Inc.

% create properties common to both patches
patchProps = {'EdgeColor','none', ...
              'FaceAlpha',0.25, ...
              'Parent',hAxes};
color = [.5 .5 .5];
yData = yLim([1 1 2 2]);

% shade the left region (if any)
if xLim(1)<xRange(1)
  xData = xscale * [xLim(1) xRange(1) xRange(1) xLim(1)];
  patch(xData, yData, color, patchProps{:});
end
  
% shade the right region (if any)
if xRange(2)<xLim(2)
  xData = xscale * [xRange(2) xLim(2) xLim(2) xRange(2)];
  patch(xData, yData, color, patchProps{:});
end