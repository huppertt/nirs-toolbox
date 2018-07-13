function show_roi_map(label,probe)

if(nargin==0)
    if(nargout==0)
        disp(nirs.util.listAtlasRegions)
    else
        varargout{1}=nirs.util.listAtlasRegions;
    end
    return
end

if(~iscellstr(label))
    label=cellstr(label);
end
label=lower(label);

aal=load(which('ROI_MNI_V5_Border_modified.mat'));
aal.BORDER_XYZ(1,:)=aal.BORDER_XYZ(1,:)*2-90;
aal.BORDER_XYZ(2,:)=aal.BORDER_XYZ(2,:)*2-126;
aal.BORDER_XYZ(3,:)=aal.BORDER_XYZ(3,:)*2-72;
aal.BORDER_XYZ=icbm_spm2tal(aal.BORDER_XYZ')';

aal.BORDER_XYZ(1,:)=aal.BORDER_XYZ(1,:)-2.5;
aal.BORDER_XYZ(2,:)=aal.BORDER_XYZ(2,:)+17.5;
aal.BORDER_XYZ(3,:)=aal.BORDER_XYZ(3,:)-20;
% 
% T =[ 0.9964    0.0178    0.0173   -0.0000
%    -0.0169    0.9957   -0.0444   -0.0000
%    -0.0151    0.0429    1.0215    0.0000
%    -0.4232  -17.5022   11.6967    1.0000];
% 
% 
% aal.BORDER_XYZ(4,:)=1;
% aal.BORDER_XYZ=(aal.BORDER_XYZ'*T)';
% aal.BORDER_XYZ=aal.BORDER_XYZ(1:3,:);
aal.ROI={aal.ROI{1} aal.ROI{2}}; 

if(nargin==2 && isa(probe,'nirs.core.Probe1020'))
    colin.mesh=probe.getmesh;
else
    colin=nirs.registration.Colin27.BEM;
end
colin.mesh(1).fiducials.Draw(:)=false;





alllabels=label;
for idx=1:length(label)
    if(isempty(strfind(label{idx},'_r')) & isempty(strfind(label{idx},'_l')))
        alllabels={alllabels{:} [label{idx} '_l'] [label{idx} '_r']};
    end
end

for i=1:length(aal.ROI)
    Labels{i}=lower(strvcat(aal.ROI{i}.Nom_L));
    Idx{i}=vertcat(aal.ROI{i}.ID);
    
    if(any(strcmp(label,'?') | strcmp(label,'*') | strcmp(label,'any')))
        lst{i}=[1:size(Labels{i},1)];
    else
        lst{i}=find(ismember(Labels{i},lower(alllabels)));
    end
end

lstNodes=[];

% Add any MNI coordinates
for i=1:length(label)
    pt=sscanf(label{i},'[%d %d %d]')';
    if(~isempty(pt))
        lstNodes(end+1)=dsearchn(aal.BORDER_XYZ',pt);
    end
end



if(isempty(lst) & isempty(lstNodes))
    disp('region not found');
    disp('use command:')
    disp(' >> nirs.util.listAtlasRegions')
    return
end


pts=colin.mesh(3).nodes;
cm=zeros(size(pts,1),1);
[k,d] = dsearchn(aal.BORDER_XYZ',pts);
cnt=1;
for i=1:length(lst)
    for j=1:length(lst{i})
        lstpt=find(aal.BORDER_V(i,k)==aal.ROI{i}(lst{i}(j)).ID);
        cm(lstpt)=cnt;
        lab{cnt}=aal.ROI{i}(lst{i}(j)).Nom_L;
        cnt=cnt+1;
        
    end
end

figure;

if(nargin>1)
    l=probe.draw3d;
    set(l,'color','k');
    hold on;
end

colin.mesh(3).transparency=1;
h=colin.mesh(3).draw(cm);
colormap(colorcube)
set(h,'FaceColor','flat')

cmap=colorcube(cnt+8);
cmap=cmap(1:cnt,:);
cmap(1,:)=.9;
colormap(cmap);
caxis([0 cnt-1])
cb=colorbar;
set(cb,'Ticks',[1:cnt-1]);
set(cb,'TickLabels',lab)
axis off;
set(gcf,'color','w')
set(cb,'TickLabelInterpreter','none')
    