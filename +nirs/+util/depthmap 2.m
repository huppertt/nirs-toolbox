function varargout = depthmap(label,headshape,atlas)

if(nargin<3)
    atlas={'aal','Brodmann (MRIcron)'};

end


if(nargin>1)
    if(isa(headshape,'nirs.core.Probe1020'))
        probe1020=headshape;
        headshape=probe1020.get_headsize;
        mesh=probe1020.getmesh;
         mesh=mesh(end);
    else
        probe1020=nirs.core.Probe1020([],headshape);
         mesh=nirs.registration.Colin27.mesh_V2;
          mesh=mesh(end);
       %  mesh=nirs.util.register_headsize
    end
else
    probe1020=nirs.core.Probe1020;
    mesh=nirs.registration.Colin27.mesh_V2;
    mesh=mesh(end);
end


if(nargin==0)
    varargout{1}=nirs.util.listAtlasRegions(mesh.labels);
    if(nargout==0)
        disp(varargout{1})
    end
    return
end


if(~iscellstr(label) && ~isnumeric(label))
    label=cellstr(label);
end

for i=1:length(label)
    if(isempty(strfind(label{i},'_R')) & isempty(strfind(label{i},'_L')))
        label{end+1}=[label{i} '_R'];
        label{end+1}=[label{i} '_L'];
    end
end

lst=find(ismember(lower(mesh.labels.keys),lower(label)));
for i=1:length(lst)
    atlas={atlas{:} mesh.labels.keys{lst(i)}};
    l=mesh.labels(mesh.labels.keys{lst(i)});
    label={label{:} l.Label{:}};
end
label={label{~ismember(label,mesh.labels.keys)}};
label=lower(label);
label=unique(label);
atlas=unique(atlas);

alllabels.Label={};
alllabels.VertexIndex={};
keys=mesh.labels.keys;
keys={keys{ismember(keys,atlas)}};

for i=1:length(keys)
    l=mesh.labels(keys{i});
    for j=1:length(l.Label)
        if(~isempty(l.Label{j}))
        alllabels.Label{end+1,1}=lower(l.Label{j});
        alllabels.VertexIndex{end+1,1}=l.VertexIndex{j};
        end
        
    end
end

lst=[];
for i=1:length(label)
  
    if(any(strcmp(label{i},'?') | strcmp(label{i},'*') | strcmp(label{i},'any')))
        lst=[lst; 1:length(alllabels.Label)];
    elseif(~isempty(strfind(label{i},'[')))
        pt=sscanf(label{i},'[%d %d %d]')';
        [k,d]=dsearchn(mesh.nodes,pt);
        alllabels.Label{end+1,1}=label{i};
        alllabels.VertexIndex{end+1,1}=k(d<5);
        lst=[lst; length(alllabels.Labels)];
    else
        lst=[lst; find(ismember(lower(alllabels.Label),label{i}))];
    end
   
end 

lst2=[];
for i=1:length(lst)
    if(isempty(alllabels.VertexIndex{lst(i)}))
        lst2=[lst2 i];
    end
end
lst(lst2)=[];


lst=unique(lst);

if(isempty(lst))
    disp('region not found');
    disp('use command:')
    disp(' >> nirs.util.listAtlasRegions')
    depth=[];
    return
end
    
headshape=probe1020.get_headsize;
tbl=nirs.util.list_1020pts('?');
Pos=[tbl.X tbl.Y tbl.Z];
  
for i=1:length(lst)
    region{i,1}=alllabels.Label{lst(i)};
    
    [~,depth(:,i)]=dsearchn(mesh.nodes(alllabels.VertexIndex{lst(i)},:),Pos);
end
[depth,i]=min(depth,[],2);
region={region{i}}';


if(nargout==0)
    figure;
    tbl2=nirs.util.list_1020pts('Cz');
    [xx,yy]=probe1020.convert2d([tbl2.X tbl2.Y tbl2.Z]);
    shiftx=-xx; shifty=-yy;
    
    [xy(:,1),xy(:,2)]=probe1020.convert2d(Pos);
    probe1020.draw1020([],[],gca);
    
    dx=mean(diff(sort(xy(:,1))));
    dy=mean(diff(sort(xy(:,2))));
    [X,Y]=meshgrid(min(xy(:,1)):dx:max(xy(:,1)),min(xy(:,2)):dy:max(xy(:,2)));
  
    warning('off','MATLAB:griddata:DuplicateDataPoints');
    IM = griddata(xy(:,1),xy(:,2),depth,X,Y,'cubic');
    
    h=imagesc([min(xy(:,1)):dx:max(xy(:,1))]+shiftx,[min(xy(:,2)):dy:max(xy(:,2))]+shifty,...
        IM,[0 max(abs(IM(:)))]);
    set(h,'alphaData',1*(~isnan(IM)));
    
    hold on;
    l=probe1020.draw1020([],[],gca);
    set(l,'LineStyle', '-', 'LineWidth', 2)
    
    set(gcf,'color','w');
    
    
    %
    %
    
    %
    if(ismember('10-20',label) | ismember('10-10',label))
        if(ismember('10-20',label))
            lst=find(ismember(tbl.Type,'10-20'));
        else
            lst=1:height(tbl);
        end
        for i=1:length(lst)
            s(i)=text(xy(lst(i),1)+shiftx,xy(lst(i),2)+shifty,tbl.Name{lst(i)});
            %         set(s(i),'Userdata',tbl.Name{i});
            %         set(s(i),'ButtonDownFcn',@displabel);
        end
        set(s,'HorizontalAlignment','center','VerticalAlignment','baseline')
    end
    
    axis tight;
    axis equal;
    axis off;
    
    cb=colorbar('SouthOutside');
    caxis([0 40])
    l=get(cb,'TickLabels');
    l{end}=['>' num2str(l{end})];
    set(cb,'TickLabels',l);
    
    
end

   

    
if(nargout>0)
    if(~isempty(probe1020.optodes_registered))
        
        %Add the link points too (since these are more useful in labels)
        ml=unique([probe1020.link.source probe1020.link.detector],'rows');
        len = size(ml, 1);
        Name = cell(len, 1); Type = Name; Units = Name;
        X = nan(len, 1); Y = X; Z = X;
        for id=1:len
            sIdx=['0000' num2str(ml(id,1))];
            sIdx=sIdx(end-3:end);
            dIdx=['0000' num2str(ml(id,2))];
            dIdx=dIdx(end-3:end);
            Name{id,1}=['Source' sIdx ':Detector' dIdx];
            X(id,1)=.5*(probe1020.swap_reg.srcPos(ml(id,1),1)+...
                probe1020.swap_reg.detPos(ml(id,2),1));
            Y(id,1)=.5*(probe1020.swap_reg.srcPos(ml(id,1),2)+...
                probe1020.swap_reg.detPos(ml(id,2),2));
            Z(id,1)=.5*(probe1020.swap_reg.srcPos(ml(id,1),3)+...
                probe1020.swap_reg.detPos(ml(id,2),3));
            
            
            %Correct for the curvature
            r=.5*(norm(probe1020.swap_reg.srcPos(ml(id,1),:))+norm(probe1020.swap_reg.detPos(ml(id,2),:)));
            theta=acos(dot(probe1020.swap_reg.srcPos(ml(id,1),:),probe1020.swap_reg.detPos(ml(id,2),:))/...
                (norm(probe1020.swap_reg.srcPos(ml(id,1),:))*norm(probe1020.swap_reg.detPos(ml(id,2),:))));
            d=r*(1-cos(theta/2));
            n=norm([X(id,1) Y(id,1) Z(id,1)]);
            X(id,1)=X(id,1)+d*X(id,1)/n;
            Y(id,1)=Y(id,1)+d*Y(id,1)/n;
            Z(id,1)=Z(id,1)+d*Z(id,1)/n;
            
            Type{id,1}='Link';
            Units{id,1}=probe1020.optodes_registered.Units{1};
        end
        probe1020.optodes_registered=[probe1020.optodes_registered;...
            table(Name,X,Y,Z,Type,Units)];
        Pts=[probe1020.optodes_registered.X ...
            probe1020.optodes_registered.Y probe1020.optodes_registered.Z];
        
        
        tbl=table;
        region={}; depth=[];
        for i=1:length(lst)
            [~,depth(:,i)]=dsearchn(mesh.nodes(alllabels.VertexIndex{lst(i)},:),Pts);
            region{i,1}=alllabels.Label{lst(i)};
        end
        [depth,i]=min(depth,[],2);
        region={region{i}}';
        depth=[probe1020.optodes_registered table(depth,region,'VariableNames',{'depth','region'})];
    else
        
        depth=[tbl table(depth,region,'VariableNames',{'depth','region'})];
    end
    varargout{1}=depth;
    
    return;
end
return


function displabel(varargin)
        
        legend(get(varargin{1},'Userdata'))
        
return



