function position = dropdownMenuAutoWidth(h)
%dropdownMenuAutoWidth sets the best width of a 'popupmenu' style uicontrol.
%  dropdownMenuAutoWidth(handle) automatically sets the best width of a
%  uicontrol of 'popupmenu' style so that every item in the drop-down list
%  is completely displayed (not being truncated).
%
%  POSITION = dropdownMenuAutoWidth(handle) returns the updated position
%  with the best width saved in POSITION(3).

%   Copyright 2003-2013 The MathWorks, Inc.

s = get(h,'String');
originalString = s;
position = get(h,'Position'); 

extents = zeros(length(originalString),4);

try
    for i = 1:length(originalString)
        extents(i,:) = get(h,'Extent');
        s = circshift(s,-1);
        set(h,'String',s);
    end
    maxWidth = max(extents(:,3));
    position(3) = maxWidth;
    set(h,'String',originalString);
catch exception
    set(h,'String',originalString);
    rethrow(exception);
end

set(h,'Position',position);