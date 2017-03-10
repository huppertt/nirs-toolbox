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
    lst=find(~any(isnan(d),1));
    [u,s,v]=nirs.math.mysvd(d(:,lst));
    proj=zeros(size(d,2),size(s,1));
    proj(lst,:)=v;
    proj=sparse(proj);
    
    data(iFile).data=u*s;
    data(iFile).projectors=proj;
    data(iFile).cov=speye(size(s,1),size(s,1));
    
    vertex=[1:size(d,2)]';
    type=repmat(cellstr('cifti'),size(d,2),1);
    
    mesh=nirs.core.Mesh;
    
    data(iFile).mesh=dtseries.core.Mesh(mesh,table(vertex,type));
end
