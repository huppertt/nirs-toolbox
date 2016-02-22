function lightsuface
% This function places a camlight in the current figure;

ax=findobj('type','axes','parent',gcf);
for i=1:length(ax)
    l=findobj('type','light','parent',ax(i));
    if(~isempty(l))
     camlight(l,'headlight');
    else
        camlight headlight;
    end
end
