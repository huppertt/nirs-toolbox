function h = plot2d(data,adderr)
% This function plots the time-series data in a spatial arrangement
% according to the probe layout


[tform,edges]=mapprobe2edges(data.probe,[0.1 0.1 .9 .9],'2D');

utypes = unique(data.probe.link.type);

if(~iscellstr(utypes))
    for i=1:length(utypes)
        utypes2{i,1}=num2str(utypes(i));
    end
    utypes=utypes2;
end

if(nargin<2)
    adderr=false;
end

types={};
for i=1:length(utypes)
    types{i}=utypes{i}(min(strfind(utypes{i},'_'))+1:end);
end
[~,~,figIdx]=unique(types);
for idx=1:length(unique(figIdx))
    figs(idx)=figure;
    set(figs(idx),'color','w');
    aa(idx)=axes('parent',figs(idx),'units','normalized','position',[0 0 1 1]);
    if(isa(data.probe,'nirs.core.Probe1020'))
        data.probe.defaultdrawfcn='2D';
    end
    l=data.probe.draw;
    axis off
    aa(idx)=get(l(1),'parent');
    
    axis tight;
    axis(aa(idx),'off')
    set(l,'color',[.85 .85 .85]);
    set(gca,'Position',[.05 .05 .9 .9]);
   % title(utypes{idx});
    sd=get(gcf,'Position');
    sd=sd(3)/20;
    set(aa(idx),'units','pixels');
end

a=[];
for idx=1:length(utypes)
    lst=find(ismember(data.probe.link.type,utypes{idx}));
    
    figure(figs(figIdx(idx)));
    hold on;
    if(iscellstr(data.probe.link.type))
        lst=find(ismember(data.probe.link.type,utypes{idx}));
    else
        lst=find(ismember(data.probe.link.type,str2num(utypes{idx})));
    end
    for j=1:length(lst)
        dIdx=data.probe.link.detector(lst(j));
        sIdx=data.probe.link.source(lst(j));
        
        p = data.probe.srcPos(sIdx,:)-.5*(data.probe.srcPos(sIdx,:)-data.probe.detPos(dIdx,:));
        
        xl=get(aa(figIdx(idx)),'XLim');
        ip=get(aa(figIdx(idx)),'Position');
        x=(p(1)-xl(1))/(xl(2)-xl(1))*ip(3)+ip(1);
        
        yl=get(aa(figIdx(idx)),'yLim');
        ip=get(aa(figIdx(idx)),'Position');
        y=(p(2)-mean(yl))/diff(yl)*ip(4)/2+ip(4)/2+ip(2);
        
        a(j)=axes;
        set(a(j),'units','pixels','position',[x-sd/2 y-sd/2 sd sd]);
        set(a(j),'units','normalized');
        set(a(j),'Ylim',[-1E-9 1e-9]);
    end
    
    for j=1:length(lst)
        hold(a(j),'on');
        s=ones(size(data.data,1),1)*real(data.data(1,lst(j)));
        if(~isreal(data.data(:,lst(j))) & adderr)
            h{idx}(lst(j))=errorbar(a(j),data.time,real(data.data(:,lst(j)))-s,...
                imag(data.data(:,lst(j))));
        else
            h{idx}(lst(j))=plot(a(j),data.time,real(data.data(:,lst(j)))-s);
        end
        set(a(j),'xlim',[min(data.time) max(data.time)]);
        set(a(j),'ylim',[min(real(data.data(:))) max(real(data.data(:))) ]);
        axis(a(j),'off');
    end
    %title(types{idx});
    set(gcf,'Name',types{idx});
    set(gcf,'NumberTitle','off');  
end
set(aa,'units','normalized');
for i=1:length(h)
    lst=[];
    for j=1:length(h{i})
        if(strcmp(class(h{i}(j)),'matlab.graphics.GraphicsPlaceholder'))
            lst=[lst j];
        end
        
    end
    h{i}(lst)=[];
end

cm=lines(length(utypes));

for i=1:length(h)
    if(~isempty(strfind(utypes{i},'hbo')))
        set(h{i},'Color','r','LineWidth',5);
    elseif(~isempty(strfind(utypes{i},'hbr')))
        set(h{i},'Color','b','LineWidth',5);
     elseif(~isempty(strfind(utypes{i},'hbt')))
        set(h{i},'Color','g','LineWidth',5);
    else
        set(h{i},'Color',cm(i,:));
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
