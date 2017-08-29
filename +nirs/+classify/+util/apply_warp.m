function datamoved = apply_warp(data,tform)
%This function applys a source-space rotation to data

if(length(data)>1)
    datamoved=data;
    for i=1:length(data)
        datamoved(i)=nirs.classify.util.apply_warp(data(i),tform);
    end
    return
end

link=data.probe.link;
link=link(ismember(link.type,link.type(1)),:);
probe=data.probe;


minX = min(probe.optodes.X);
maxX = max(probe.optodes.X);
dX = (maxX-minX)/3;
minY = min(probe.optodes.Y);
maxY = max(probe.optodes.Y);
dY = (maxY-minY)/3;

[X,Y,Z]=meshgrid([minX-dX:dX/30:maxX+dX],[minY-dY:dY/30:maxY+dY],[-10]);

mesh=nirs.core.Mesh;
mesh.nodes=[X(:) Y(:) Z(:)];

lambda=unique(probe.link.type);
probe.link=probe.link(ismember(probe.link.type,probe.link.type(1)),:);
if(~isnumeric(lambda))
    probe.link.type=ones(height(probe.link),1)*808;
    lambda=808;
end

FwdModel=nirs.forward.ApproxSlab;

FwdModel.prop=nirs.media.tissues.brain(.7,50,lambda);
FwdModel.Fm=0;
FwdModel.probe=probe;
FwdModel.mesh=mesh;
J=FwdModel.jacobian;
J=J.mua;

for j=1:2
    moving=reshape(mesh.nodes(:,j),size(X));
    mesh.nodes(:,j) = reshape(imwarp(moving,tform.invert,'OutputView',imref2d(size(moving))),[],1);
end
FwdModel.mesh=mesh;
J2=FwdModel.jacobian;
J2=J2.mua;
a=J*pinv(J2);




datamoved=data.sorted();
if(isa(data,'nirs.core.Data') | isa(data,'eeg.core.Data'))
    n=size(datamoved.data,2)/size(a,1);
    a=kron(a,eye(n));
    datamoved.data=datamoved.data*a';
   
elseif(isa(data,'nirs.core.ChannelStats') | isa(data,'eeg.core.ChannelStats'))
    n=length(datamoved.beta)/size(a,1);
    a=kron(a,eye(n));
    datamoved.beta=a*datamoved.beta;
    datamoved.covb=a*datamoved.covb*a';
else
    error('type not supported');
end