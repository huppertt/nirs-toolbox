function lightsuface(h)
% This function places a camlight in the current figure;

if(nargin==0)
    h=gcf;
end

ax=findobj('type','axes','parent',h);
for i=1:length(ax)
    l=findobj('type','light','parent',ax(i));
    if(~isempty(l))
        if(length(l)>1)
            delete(l(2:end));
            l=findobj('type','light','parent',ax(i));
        end
     camlight(l,'headlight');
    else
        camlight headlight;
    end
end
