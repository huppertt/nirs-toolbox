function fwdModel = read_freesurfer(SubjectsDIR,lambda)


if(nargin<2)
    lambda=[690 830];
end

% this assumes MNI watershed has been run

f=[rdir(fullfile(SubjectsDIR,'bem','watershed','**','*_outer_skin_surface')); ...
    rdir(fullfile(SubjectsDIR,'bem','outer_skin.surf'))];
    
[nodes,faces]=read_surf(f(1).name);
mesh(1,1)=nirs.core.Mesh(nodes,faces+1);
mesh(1).transparency=.1;
f=[rdir(fullfile(SubjectsDIR,'bem','watershed','**','*_outer_skull_surface')) ...
    rdir(fullfile(SubjectsDIR,'bem','outer_skull.surf'))];
[nodes,faces]=read_surf(f(1).name);
mesh(2,1)=nirs.core.Mesh(nodes,faces+1);
mesh(2).transparency=.2;
f=[rdir(fullfile(SubjectsDIR,'bem','watershed','**','*_inner_skull_surface'))...
    rdir(fullfile(SubjectsDIR,'bem','inner_skull.surf'))];
[nodes,faces]=read_surf(f(1).name);
mesh(3,1)=nirs.core.Mesh(nodes,faces+1);
mesh(3).transparency=.2;

% Use the wavelet methods defined in 
% Phys Med Biol. 2009 Oct 21;54(20):6383-413. doi: 10.1088/0031-9155/54/20/023. Epub 2009 Oct 7.
% Topographic localization of brain activation in diffuse optical imaging using spherical wavelets.
% Abdelnour F1, Schmidt B, Huppert TJ.
J=6;  %levels of wavelets in the model (=ico6)

options=struct('base_mesh', 'ico', 'keep_subdivision', true);
[vertex,face] = compute_semiregular_sphere(J,options);
 
[lsv,lsf]=read_surf(fullfile(SubjectsDIR,'surf','lh.sphere.reg'));
[rsv,rsf]=read_surf(fullfile(SubjectsDIR,'surf','rh.sphere.reg'));

kl=dsearchn(lsv/100,vertex{end}');
kr=dsearchn(rsv/100,vertex{end}');

% Create the surfaces
[lv,lf]=read_surf(fullfile(SubjectsDIR,'surf','lh.pial'));
[rv,rf]=read_surf(fullfile(SubjectsDIR,'surf','rh.pial'));

mesh(4,1)=nirs.core.Mesh([lv(kl,:); rv(kr,:)],[face{end}'; face{end}'+length(kl)]);  



% get the fiducials



%%  Register the 10-5 labels to the head
% These are the values in 
% /home/pkg/software/MNE/share/mne/mne_analyze/fsaverage/fsaverage-fiducials.fif
% It just doesn't make sense have to read each time so I hard coded it
fid =[1.5000   85.1000  -34.8000    1.0000   % nas
    -80.6000  -29.1000  -41.3000    1.0000   % lpa
     84.4000  -28.5000  -41.3000    1.0000]; %rpa

T=importdata(fullfile(SubjectsDIR,'mri','transforms','talairach.xfm'));
T.data(4,4)=1;
T=T.data;
pts=fid*T';
pts(:,4)=[];

% Add the 10-20 fiducials to the scalp mesh

tbl1020=nirs.util.list_1020pts;
xyz1020=[tbl1020.X tbl1020.Y tbl1020.Z];



a=pts(1,:)-.5*(pts(2,:)+pts(3,:));
b=pts(2,:)-.5*(pts(2,:)+pts(3,:));
a=a/norm(a);
b=b/norm(b);
v=cross(a,b);
v=v/norm(v);
x=ones(4001,1)*(.5*(pts(2,:)+pts(3,:)))+[0:.05:200]'*v;

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
x=ones(4001,1)*(.5*(pts2(2,:)+pts2(3,:)))+[0:.05:200]'*v;
[k,d]=dsearchn(xyz1020,x);
[~,id]=min(d);
pts2(4,:)=x(id,:);

pts(:,4)=1;
pts2(:,4)=1;
xyz1020(:,4)=1;

T=pts2\pts;
xyz1020=xyz1020*T;  
xyz1020(:,4)=[];

% now registered to the head
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

fi={'aparc.a2009s.annot','aparc.DKTatlas40.annot'};

Labels=Dictionary;
    % add the Atlas labels
for id=1:length(fi)
    [Rvertices, Rlabel, Rcolortable] =read_annotation(fullfile(SubjectsDIR,'label',['rh.' fi{id}]));
    [Lvertices, Llabel, Lcolortable] =read_annotation(fullfile(SubjectsDIR,'label',['lh.' fi{id}]));
    
    [~,a]=ismember(kr,Rvertices+1);
    RlabelsIdx=Rlabel(a);
    [~,a]=ismember(kl,Lvertices+1);
    LlabelsIdx=Llabel(a);
    
    
    
    
    S=struct('VertexIndex',{''},'Label',{''},'Region',{''});
    cnt=1;
    for i=1:length(Lcolortable.struct_names)
        lst=find(LlabelsIdx==Lcolortable.table(i,5));
        if(length(lst)>0)
            name=['L_' Lcolortable.struct_names{i}];
            S.VertexIndex{cnt,1}=lst;
            S.Label{cnt,1}=name;
            S.Region{cnt,1}='left cortex';
            cnt=cnt+1;
        end
    end
    for i=1:length(Rcolortable.struct_names)
        lst=find(RlabelsIdx==Rcolortable.table(i,5));
        if(length(lst)>0)
            name=['R_' Lcolortable.struct_names{i}];
            S.VertexIndex{cnt,1}=lst+length(kl);
            S.Label{cnt,1}=name;
            S.Region{cnt,1}='right cortex';
            cnt=cnt+1;
        end
    end
    Labels(fi{id})=S;
end


    
mesh(4).labels=Labels;

% Now make the forward model
fwdModel=nirs.forward.NirfastBEM;
fwdModel.mesh=mesh;
fwdModel.prop={nirs.media.tissues.skin(lambda)...
               nirs.media.tissues.bone(lambda)...
               nirs.media.tissues.water(lambda)...
               nirs.media.tissues.brain(lambda,.7,60)};