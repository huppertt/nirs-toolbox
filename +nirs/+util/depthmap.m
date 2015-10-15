function depth = depthmap(label,Pos);

if(~iscellstr(label))
    label=cellstr(label);
end
label=lower(label);

aal=load(which('ROI_MNI_V5_Border.mat'));
aalLabels=load(which('ROI_MNI_V5_List.mat'));
aal.BORDER_XYZ(1,:)=aal.BORDER_XYZ(1,:)*2-90;
aal.BORDER_XYZ(2,:)=aal.BORDER_XYZ(2,:)*2-126;
aal.BORDER_XYZ(3,:)=aal.BORDER_XYZ(3,:)*2-72;

if(nargin<2 || isempty(Pos))
    fid=fopen(which('ext1020.sfp'),'r');
    marker=textscan(fid,'%s\t%d\t%d\t%d');
    fclose(fid);
    Pos=double([marker{2} marker{3} marker{4}]);
end

alllabels=label;
for idx=1:length(label)
    if(isempty(strfind(label{idx},'_r')) & isempty(strfind(label{idx},'_l')))
        alllabels={alllabels{:} [label{idx} '_l'] [label{idx} '_r']};
    end
end


Labels=lower(strvcat(aalLabels.ROI.Nom_L));
Idx=vertcat(aalLabels.ROI.ID);
lst=find(ismember(Labels,lower(alllabels)));

if(isempty(lst))
    disp('region not found');
    disp('use command:')
    disp(' >> nirs.util.listAtlasRegions')
    depth=[];
    return
end

lstNodes=find(ismember(aal.BORDER_V,Idx(lst)));

[k,depth] = dsearchn(aal.BORDER_XYZ(:,lstNodes)',Pos);

if(nargout==0)
    figure;
    headradius = mean(sqrt(sum(Pos.^2,2)));
    
    %Clarke far-side general prospective azumuthal projection
    r = -2.4;
    R=sqrt(sum(Pos.^2,2));
    xy(:,1)=r*R.*(Pos(:,1)./abs(Pos(:,3)-r*R));
    xy(:,2)=r*R.*(Pos(:,2)./abs(Pos(:,3)-r*R));
    
    
    dx=mean(diff(sort(xy(:,1))));
    dy=mean(diff(sort(xy(:,2))));
    [X,Y]=meshgrid(min(xy(:,1)):dx:max(xy(:,1)),min(xy(:,2)):dy:max(xy(:,2)));
    warning('off','MATLAB:griddata:DuplicateDataPoints');
    IM = griddata(xy(:,1),xy(:,2),depth,X,Y,'cubic');
    
    h=imagesc(min(xy(:,1)):dx:max(xy(:,1)),min(xy(:,2)):dy:max(xy(:,2)),...
        IM,[0 max(abs(IM(:)))]);
    set(h,'alphaData',1*(~isnan(IM)));
    
    hold on;
    set(gcf,'color','w');
    
    scatter(xy(:,1),xy(:,2),'filled','MarkerFaceColor',[.8 .8 .8]);
    
    axis tight;
    axis equal;
    axis off;
    
    % add a circle for the head
    theta = linspace(0,2*pi);
    plot(headradius*cos(theta),headradius*sin(theta),'k');
    line([-10 0],[-headradius -headradius-10],'color','k');
    line([10 0],[-headradius -headradius-10],'color','k');
    scatter([-15 15],[-headradius+5 -headradius+5],'filled','k','sizedata',120);
    set(gca,'YDir','reverse');
    colorbar;
    caxis([0 30])
end