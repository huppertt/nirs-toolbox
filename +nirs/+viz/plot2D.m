function h = plot2d(data,samefig)
% This function plots the time-series data in a spatial arrangement
% according to the probe layout


[tform,edges]=mapprobe2edges(data.probe,[0.1 0.1 .9 .9],'2D');

utypes = unique(data.probe.link.type);

if(nargin<2)
    samefig=false;
end

a=[];
for idx=1:length(utypes)
    lst=find(ismember(data.probe.link.type,utypes{idx}));
    if(~samefig || isempty(a))
        f=figure;
        set(f,'color','w');
        
        
        aa=axes('parent',f,'units','normalized','position',[0 0 1 1]);
        l=data.probe.draw;
        axis tight;
        axis(aa,'off')
        set(l,'color',[.85 .85 .85]);
        set(gca,'Position',[0.1 0.1 .8 .8]);
        title(utypes{idx});
        
        for j=1:length(lst)
            dIdx=data.probe.link.detector(lst(j));
            sIdx=data.probe.link.source(lst(j));
            
            p = .5*(data.probe.srcPos(sIdx,:)+data.probe.detPos(dIdx,:));
            p(:,3)=1;
            p=p*tform;
            pos=[p(1)-edges(1) p(2)-edges(2) edges(1)*2 edges(2)*2];
            a(j)=axes('parent',f,'units','normalized','position',pos);
        end
        lstA=lst;
        
      
    else
        lstA=[1:size(data.data,2)];
        
        title(aa,strvcat(utypes{:}));
    end
    for j=1:length(lst)
        hold(a(j),'on');
        plot(a(j),data.time,data.data(:,lst(j)));
        set(a(j),'xlim',[min(data.time) max(data.time)]);
        set(a(j),'ylim',[min(min(data.data(:,lstA))) max(max(data.data(:,lstA)))]);
        axis(a(j),'off');
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
