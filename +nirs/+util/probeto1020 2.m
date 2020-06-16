function probe1020=probeto1020(probe)

probe1020=nirs.core.Probe1020;



% Get the 10-20 labels from the SPM template
fid=fopen(which('ext1020.sfp'),'r');
marker=textscan(fid,'%s\t%d\t%d\t%d');
fclose(fid);
Pos1020=double([marker{2} marker{3} marker{4}]);

headradius = mean(sqrt(sum(Pos1020.^2,2)));

%Clarke far-side general prospective azumuthal projection 
r = -2.4;
R=sqrt(sum(Pos1020.^2,2));
xy1020(:,1)=r*R.*(Pos1020(:,1)./abs(Pos1020(:,3)-r*R));
xy1020(:,2)=r*R.*(Pos1020(:,2)./abs(Pos1020(:,3)-r*R));

% Points from the probe
Pos(:,1)=probe.optodes.X;
Pos(:,2)=probe.optodes.Y;
Pos(:,3)=probe.optodes.Z;
R=sqrt(sum(Pos.^2,2));
xy(:,1)=r*R.*(Pos(:,1)./abs(Pos(:,3)-r*R));
xy(:,2)=r*R.*(Pos(:,2)./abs(Pos(:,3)-r*R));

if(~iscellstr(probe.link.type))
    probe.link.type=arrayfun(@(x)(cellstr(num2str(x))),probe.link.type);
end

utype=unique(probe.link.type);

figure;
for u=1:length(utype)
    subplot(1,length(utype),u);
    
    lst=find(ismember(probe.link.type,utype{u}));
    
    if(~isempty(im))
        warning('need to fix this');
        dx=mean(diff(sort(xy(lst,1))));
        dy=mean(diff(sort(xy(lst,2))));
        [X,Y]=meshgrid(min(xy(lst,1)):dx:max(xy(lst,1)),min(xy(lst,2)):dy:max(xy(lst,2)));
        IM = griddata(xy(lst,1),xy(lst,2),im(lst),X,Y,'cubic');
        
        h=imagesc(min(xy(lst,1)):dx:max(xy(lst,1)),min(xy(lst,2)):dy:max(xy(lst,2)),...
            IM,[-max(abs(IM(:))) max(abs(IM(:)))])
        set(h,'alphaData',1*(~isnan(IM)));
    end
    
    hold on;
    set(gcf,'color','w')
    
    scatter(xy1020(:,1),xy1020(:,2),'filled','MarkerFaceColor',[.8 .8 .8])
    
    lstS=find(ismember(probe.optodes.Type,'Source'));
    scatter(xy(lstS,1),xy(lstS,2),'filled','MarkerFaceColor','r')
    lstD=find(ismember(probe.optodes.Type,'Detector'));
    scatter(xy(lstD,1),xy(lstD,2),'filled','MarkerFaceColor','b')
    
    for i=1:height(probe.link)
        s=probe.link.source(i);
        d=probe.link.detector(i);
        l(i)=line(xy([lstS(s) lstD(d)],1),xy([lstS(s) lstD(d)],2),'color','k');
    end
    
    axis tight
    axis equal
    axis off;
    
    % add a circle for the head
    theta = linspace(0,2*pi);
    plot(headradius*cos(theta),headradius*sin(theta),'k');
    line([-10 0],[-headradius -headradius-10],'color','k')
    line([10 0],[-headradius -headradius-10],'color','k')
    scatter([-15 15],[-headradius+5 -headradius+5],'filled','k','sizedata',120)
    
    set(gca,'YDir','reverse')
    
end