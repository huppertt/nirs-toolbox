function data = loadCifti(filenames)

% if a single filename, put it in a cell
if ischar( filenames )
    filenames = {filenames};
end

data = dtseries.core.Data;
data(:)=[];

% iterate through cell array
for iFile = 1:length(filenames)
    c=ft_read_cifti(filenames{iFile});
    data(iFile).description=filenames{iFile};
    data(iFile).time=c.time;
    
    lst=find(ismember(c.brainstructure,find(ismember(c.brainstructurelabel,{'CORTEX_LEFT','CORTEX_RIGHT'}))));
    d=c.dtseries(lst,:)';
    
    vertex=[1:size(d,2)]';
    lst=find([any(isnan(d),1) | all(d==0,1) | sqrt(var(d,[],1))<eps(1)*10]);
    vertex(lst)=[];
    d(:,lst)=[];
    
    [u,s,proj]=nirs.math.mysvd(d);
    lst=find(diag(s)<eps(single(1)));
    u(:,lst)=[]; s(lst,:)=[]; s(:,lst)=[]; v(lst,:)=[];
    
    proj=sparse(proj);
    
    data(iFile).data=u*s;
    data(iFile).projectors=proj;
    data(iFile).cov=speye(size(s,1),size(s,1));
    
    
    type=repmat(cellstr('cifti'),size(d,2),1);
    
    mesh=nirs.core.Mesh;
    
    data(iFile).mesh=dtseries.core.Mesh(mesh,table(vertex,type));
end
