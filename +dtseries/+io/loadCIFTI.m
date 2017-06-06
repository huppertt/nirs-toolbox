function data = loadCifti(filenames,surf)
% p='/Users/huppert/Desktop/Cerebro/1_R21_Data/2_fNIRS-EEG-fMRI/fMRI/analyzed/NVC002/MNINonLinear/Results/ep2d_bold_MN1';
% p2='/Users/huppert/Desktop/Cerebro/1_R21_Data/2_fNIRS-EEG-fMRI/fMRI/analyzed/NVC002/MNINonLinear/fsaverage_LR32k';
% filenames={fullfile(p,'ep2d_bold_MN1_Atlas.dtseries.nii')};
% surf = {fullfile(p2,'NVC002.L.pial.32k_fs_LR.surf.gii')};


% if a single filename, put it in a cell
if ischar( filenames )
    filenames = {filenames};
end
if ischar( surf )
   surf = {surf};
end


data = dtseries.core.Data;
data(:)=[];

% iterate through cell array
for iFile = 1:length(filenames)
    c=ft_read_cifti(filenames{iFile});
    data(iFile).description=filenames{iFile};
    data(iFile).time=c.time;
    
    lst=find(ismember(c.brainstructure,find(ismember(c.brainstructurelabel,{'CORTEX_LEFT','CORTEX_RIGHT'}))));
    
    [p,surfroot,ext]=fileparts(surf{iFile});
    surfL=surfroot;
    surfL([strfind(surfL,'.L.') strfind(surfL,'.R.')]+1)='L';
    surfR=surfroot;
    surfR([strfind(surfR,'.L.') strfind(surfR,'.R.')]+1)='R';
    giiL=gifti(fullfile(p,[surfL ext]));
    giiR=gifti(fullfile(p,[surfR ext]));
     
    v = [giiL.vertices; giiR.vertices];
    f = [giiL.faces; giiR.faces+size(giiL.vertices,1)];
    
    d=c.dtseries(lst,:)';
    
    vertex=[1:size(d,2)]';
    lst=find([any(isnan(d),1) | all(d==0,1) | sqrt(var(d,[],1))<eps(1)*10]);
    vertex(lst)=[];
    d(:,lst)=[];
    
    [u,s,proj]=nirs.math.mysvd(d);
    lst=find(diag(s)<eps(single(1)));
    u(:,lst)=[]; s(lst,:)=[]; s(:,lst)=[]; proj(lst,:)=[];
    
    proj=sparse(proj);
    
    data(iFile).data=u*s;
    data(iFile).projectors=proj;
    data(iFile).cov=speye(size(proj,2),size(proj,2));
    
    
    type=repmat(cellstr('cifti'),size(d,2),1);
    
    mesh=nirs.core.Mesh(v,f);
    
    data(iFile).mesh=dtseries.core.Mesh(mesh,table(vertex,type));
end
