function plotpkmarkers(hLine, y)
%PLOTPKMARKERS  handler to offset peak markers by the size of the marker
%
%   This function is for internal use only.
%
%   If the figure zoom ActionPostCallback event cannot be overridden
%   warn user that we can't override and use (centered) circular markers
%   for peaks.
%
%   Otherwise, use inverted triangular markers and offset them by the size
%   of the marker, taking care to recompute the offsets upon size change
%   and figure zoom.

%   Copyright 2014 The MathWorks, Inc.

if ~isempty(y) && ~isempty(hLine) && ishghandle(hLine)
  hFig = ancestor(hLine,'figure');
  hZoom = zoom(hFig);  
  cbActionPostZoom = get(hZoom,'ActionPostCallback');
  cbZoom = createZoomCallback();
  if isempty(cbActionPostZoom) || isequal(char(cbActionPostZoom),char(cbZoom))
    % offset markers when figure zoom is invoked.
    if isempty(cbActionPostZoom)
      set(hZoom, 'ActionPostCallback', cbZoom);
    end
    
    % use solid inverted triangular markers
    set(hLine, 'Marker', 'v');
    set(hLine, 'MarkerFaceColor', get(hLine,'Color'));

    % create callback to offset the peak markers
    setappdata(hLine, 'YPos', y);
    cbSizeChange = createSizeChangeCallback(hLine);
    
    % offset markers on a SizeChange/SizeChanged event.
    listener = event.listener( hFig, 'SizeChanged', cbSizeChange);
    setappdata(hLine, 'Listener', listener);

    %Offset the markers now.
    cbSizeChange(hLine);
  else
    % zoom post-action callback has already been overridden
    warning(message('signal:findpeaks:CantOffsetPeakMarkers'));
  end
end

function fcn = createSizeChangeCallback(hLine)
fcn = @(obj, evt) offsetPeakMarkers(hLine);

function fcn = createZoomCallback()
fcn = @(hFig, hAxes) offsetPeakMarkersByAxes(hFig, hAxes);

% offset all peak markers contained within an axes
function offsetPeakMarkersByAxes(~, evt)
if isfield(evt,'Axes') && ishghandle(evt.Axes)
  hLines = findall(evt.Axes,'Type','line','Tag','Peak');
  for i=1:numel(hLines)
    offsetPeakMarkers(hLines(i));
  end
end

% offset the markers in the line
function offsetPeakMarkers(hLine)
if ishghandle(hLine)
  % Fetch the data needed to compute the offset line
  y = getappdata(hLine, 'YPos');
  hAxes = ancestor(hLine, 'axes');
  yMarkerOffset = get(hLine,'MarkerSize');
  axesPos = getpixelposition(hAxes);
  yLim = get(hAxes,'YLim');

  % bump the line y data by the marker size
  yOffset = yMarkerOffset * diff(yLim) ./ axesPos(4);
  yOld = get(hLine, 'YData');
  yNew = y + yOffset;
  if ~isequaln(yOld, yNew)
    set(hLine,'YData', y + yOffset);
  end
end
