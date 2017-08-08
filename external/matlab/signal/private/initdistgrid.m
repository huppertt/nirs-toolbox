function initdistgrid(hAxes)
%INITDISTGRID initialize distortion plot grid
%   Initialize the distortion grid so that it uses a light gray grid.
%   Additionally make room for markers.
%   
%   This function is for internal use only. It may be removed in the future.
   
%   Copyright 2013-2014 The MathWorks, Inc.

% bound plot by y limit; give 10 extra dB for fundamental marker.
yTick = get(hAxes,'YTick');
if numel(yTick>1)
  yNew = yTick(end) + 10;
  yLim = get(hAxes,'YLim');
  set(hAxes,'YLim',[yLim(1) yNew]);
end

% Ensure axes limits are properly cached for zoom/unzoom
resetplotview(hAxes,'SaveCurrentView');  

% turn on box and grid
set(hAxes,'Box','on','XGrid','on','YGrid','on');