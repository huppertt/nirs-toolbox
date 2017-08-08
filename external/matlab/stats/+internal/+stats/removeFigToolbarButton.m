function removeFigToolbarButton(figureHandle,removeTags)
%REMOVEFIGTOOLBARBUTTON Removes unnecessary buttons from the Figure
%Toolbar. Buttons are specified by their Tags. removeTags is a cell array
%of strings.
%   Example: 
%   To Removes the Brushing, InsertColorBar and Rotate buttons.
%   >> internal.stats.removeFigToolbarButton(gcf,{'Exploration.Brushing',...
%                                           'Annotation.InsertColorbar',...
%                                            'Exploration.Rotate'});
%   

%   Copyright 2012 The MathWorks, Inc.


h0 = findall(figureHandle,'Type','uitoolbar');
h1 = findall(h0,'Parent',h0);

btnTags = get(h1,'Tag');

for i = 1:numel(btnTags)
    if any(strcmp(btnTags{i},removeTags))
        delete(h1(i));
    end
end
end