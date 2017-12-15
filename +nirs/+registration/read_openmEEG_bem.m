function fwdModel = read_openmEEG_bem(filename,lambda)
% Conversion code for reading OpenMEEG - BrainStorm files

%filename='brainstormsubject.mat';

if(nargin<2)
    lambda=[690 830];
end


f=load(filename);

fields={'Anatomy','Cortex','Scalp','OuterSkull','InnerSkull'};

pth=fileparts(filename);

for i=1:length(fields)
    if(exist(f.(fields{i}),'file')==2)
        f.(fields{i})=load(f.(fields{i}));
    elseif(exist(fullfile(pth,f.(fields{i})(max(strfind(f.(fields{i}),filesep))+1:end)),'file')==2)
        %try to read from current folder
        f.(fields{i})=load(fullfile(pth,f.(fields{i})(max(strfind(f.(fields{i}),filesep))+1:end)));
    elseif(exist(f.(fields{i})(max(strfind(f.(fields{i}),filesep))+1:end),'file')==2)
        %try to read from current folder
        f.(fields{i})=load(f.(fields{i})(max(strfind(f.(fields{i}),filesep))+1:end));     
        
    else
        warning(['unable to load ' fields{i} ' file']);
        
    end
end

% Make the Mesh models
mesh(1,1)=nirs.core.Mesh(f.Scalp.Vertices*1000,f.Scalp.Faces);
mesh(1).transparency=.1;

mesh(2,1)=nirs.core.Mesh(f.OuterSkull.Vertices*1000,f.OuterSkull.Faces);
mesh(2).transparency=.2;
mesh(3,1)=nirs.core.Mesh(f.InnerSkull.Vertices*1000,f.InnerSkull.Faces);
mesh(3).transparency=.2;

mesh(4,1)=nirs.core.Mesh(f.Cortex.Vertices*1000,f.Cortex.Faces);


% Add the 10-20 fiducials to the scalp mesh

tbl1020=nirs.util.list_1020pts;
xyz1020=[tbl1020.X tbl1020.Y tbl1020.Z];

pts(1,:)=f.Scalp.SCS.NAS;
pts(2,:)=f.Scalp.SCS.LPA;
pts(3,:)=f.Scalp.SCS.RPA;
pts=(pts-ones(3,1)*f.Scalp.SCS.Origin')*f.Scalp.SCS.R';


a=pts(1,:)-.5*(pts(2,:)+pts(3,:));
b=pts(2,:)-.5*(pts(2,:)+pts(3,:));
a=a/norm(a);
b=b/norm(b);
v=cross(a,b);
v=v/norm(v);
x=.5*(pts(2,:)+pts(3,:))+[0:.05:200]'*v;

[k,d]=dsearchn(mesh(1).nodes,x);
[~,id]=min(d);
pts(4,:)=x(id,:);

pts2=[tbl1020.X(1:3) tbl1020.Y(1:3) tbl1020.Z(1:3)]; 
a=pts2(1,:)-.5*(pts2(2,:)+pts2(3,:));
b=pts2(2,:)-.5*(pts2(2,:)+pts2(3,:));
a=a/norm(a);
b=b/norm(b);
v=cross(a,b);
v=v/norm(v);
x=.5*(pts2(2,:)+pts2(3,:))+[0:.05:200]'*v;
[k,d]=dsearchn(xyz1020,x);
[~,id]=min(d);
pts2(4,:)=x(id,:);

pts(:,4)=1;
pts2(:,4)=1;
xyz1020(:,4)=1;

T=pts2\pts;
xyz1020=xyz1020*T;  % now registered to the head
xyz1020(:,4)=[];
% let's do an itertive closest point to refine
% (icp is included in my NIRS-toolbox)

for i=1:3
    [TR, TT] = icp(mesh(1).nodes',xyz1020');
    xyz1020=(TR * xyz1020' + TT*ones(1,size(xyz1020,1)))';
    k=dsearchn(xyz1020(:,1:3),mesh(1).nodes);
    xyz1020(:,1:3)=mesh(1).nodes(dsearchn(mesh(1).nodes,xyz1020(:,1:3)),:);
end

tbl1020.X=xyz1020(:,1);
tbl1020.Y=xyz1020(:,2);
tbl1020.Z=xyz1020(:,3);

mesh(1)=mesh(1).addfiducials(tbl1020);


Labels=Dictionary;
% add the Atlas labels
for i=1:length(f.Cortex.Atlas)
    if(length(f.Cortex.Atlas(i).Scouts)>0)
        name=f.Cortex.Atlas(i).Name;
        S=struct('VertexIndex',{''},'Label',{''},'Region',{''});
        for j=1:length(f.Cortex.Atlas(i).Scouts);
            S.VertexIndex{j}=f.Cortex.Atlas(i).Scouts(j).Vertices;
            S.Label{j}=f.Cortex.Atlas(i).Scouts(j).Label;
            S.Region{j}=f.Cortex.Atlas(i).Scouts(j).Region;
        end
        
    Labels(name)=S;
    end
end

mesh(4).labels=Labels;

% Now make the forward model
fwdModel=nirs.forward.NirfastBEM;
fwdModel.mesh=mesh;
fwdModel.prop={nirs.media.tissues.skin(lambda)...
               nirs.media.tissues.bone(lambda)...
               nirs.media.tissues.water(lambda)...
               nirs.media.tissues.brain(.7,60,lambda)};





