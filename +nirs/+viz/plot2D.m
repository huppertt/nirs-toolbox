function h = plot2d(data,samefig,adderr)
% This function plots the time-series data in a spatial arrangement
% according to the probe layout


[tform,edges]=mapprobe2edges(data.probe,[0.1 0.1 .9 .9],'2D');

utypes = unique(data.probe.link.type);

if(nargin<2)
    samefig=false;
end
if(nargin<3)
    adderr=false;
end


a=[];
for idx=1:length(utypes)
    lst=find(ismember(data.probe.link.type,utypes{idx}));
    if(~samefig || isempty(a))
        f=figure;
        set(f,'color','w');
        
        
        aa=axes('parent',f,'units','normalized','position',[0 0 1 1]);
        if(isa(data.probe,'nirs.core.Probe1020'))
             data.probe.defaultdrawfcn='2D';
        end
            l=data.probe.draw;
            axis off
        aa=get(l(1),'parent');
        
        axis tight;
        axis(aa,'off')
        set(l,'color',[.85 .85 .85]);
        set(gca,'Position',[.05 .05 .9 .9]);
        title(utypes{idx});
        sd=get(gcf,'InnerPosition');
        sd=sd(3)/20;
        set(aa,'units','pixels');
        hold on;
        for j=1:length(lst)
            dIdx=data.probe.link.detector(lst(j));
            sIdx=data.probe.link.source(lst(j));
            
            p = data.probe.srcPos(sIdx,:)-.5*(data.probe.srcPos(sIdx,:)-data.probe.detPos(dIdx,:));
            
            xl=get(aa,'XLim');
            ip=get(aa,'Position');
            x=(p(1)-xl(1))/(xl(2)-xl(1))*ip(3)+ip(1);
            
            yl=get(aa,'yLim');
            ip=get(aa,'Position');
            y=(p(2)-mean(yl))/diff(yl)*ip(4)/2+ip(4)/2+ip(2);
            
             a(j)=axes;
            set(a(j),'units','pixels','position',[x-sd/2 y-sd/2 sd sd]);
            set(a(j),'units','normalized');
        end
       set(aa,'units','normalized');
        lstA=lst;
        
      
    else
        lstA=[1:size(data.data,2)];
        
        %title(aa,strvcat(utypes{:}));
    end
    for j=1:length(lst)
        hold(a(j),'on');
        if(~isreal(data.data(:,lst(j))) & adderr)
            h{idx}(lst(j))=errorbar(a(j),data.time,real(data.data(:,lst(j))),...
                imag(data.data(:,lst(j))));
        else
            h{idx}(lst(j))=plot(a(j),data.time,real(data.data(:,lst(j))));
        end
        set(a(j),'xlim',[min(data.time) max(data.time)]);
        set(a(j),'ylim',[min(min(real(data.data(:,lstA)))) max(max(real(data.data(:,lstA))))]);
        axis(a(j),'off');
    end
   
end

for i=1:length(h)
    lst=[];
    for j=1:length(h{i})
        if(strcmp(class(h{i}(j)),'matlab.graphics.GraphicsPlaceholder'))
            lst=[lst j];
        end
        
    end
    h{i}(lst)=[];
end

for i=1:length(h)
    if(~isempty(strfind(utypes{i},'hbo')))
        set(h{i},'Color','r');
    elseif(~isempty(strfind(utypes{i},'hbr')))
        set(h{i},'Color','b');
     elseif(~isempty(strfind(utypes{i},'hbt')))
        set(h{i},'Color','g');
    end
end


end

function [tform,edges]=mapprobe2edges(probe,position,type)
% This function maps the edges of the probe to match the dimensions of the
% figure.
% x/y are the centers of each window
% dx/dy are 1/2 the window edges

minXY = min([probe.srcPos; probe.detPos],[],1);
maxXY = max([probe.srcPos; probe.detPos],[],1);

p(1,:)=[minXY(1) minXY(2) 1];
p(2,:)=[maxXY(1) minXY(2) 1];
p(3,:)=[maxXY(1) maxXY(2) 1];
p(4,:)=[minXY(1) maxXY(2) 1];

p2(1,:)=[position(1) position(2) 1];
p2(2,:)=[position(1)+position(3) position(2) 1];
p2(3,:)=[position(1)+position(3) position(2)+position(4) 1];
p2(4,:)=[position(1) position(2)+position(4) 1];

% Do it this way to make it easier to map to non-linear spaces later (e.g.
% 10-20 spherical plots)

tform=p\p2;

nx=length(unique([probe.srcPos(:,1); probe.detPos(:,1)]));
ny=length(unique([probe.srcPos(:,2); probe.detPos(:,2)]));

edges(1)=position(3)/(nx+1)/2;
edges(2)=position(4)/(ny+1)/2;



end
