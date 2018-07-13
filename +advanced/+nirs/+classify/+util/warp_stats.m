function [datamoved,tform]=warp_stats(datatarget,datamovable,type)

if(nargin<3)
    type='affine';
%        TRANSFORMTYPE         TYPES OF DISTORTION
%        -------------         -----------------------
%        'translation'         Translation
%        'rigid'               Translation, Rotation
%        'similarity'          Translation, Rotation, Scale
%        'affine'              Translation, Rotation, Scale, Shear
end

varmov=datamovable.variables;
vartarg=datatarget.variables;

varmov=varmov(:,1:2);
vartarg=vartarg(:,1:2);

[i1,~,j1]=unique(varmov);
[i2,~,j2]=unique(vartarg);

datatarg=zeros(height(i1),length(find(j1==1)));
datamov=zeros(height(i2),length(find(j2==1)));
datatargT=zeros(height(i1),length(find(j1==1)));
datamovT=zeros(height(i2),length(find(j2==1)));

for i=1:height(i1)
    datatarg(i,:)=datatarget.beta(j1==i);
    datatargT(i,:)=datatarget.tstat(j1==i);
end

for i=1:height(i2)
    datamov(i,:)=datamovable.beta(j2==i);
    datamovT(i,:)=datamovable.tstat(j2==i);
end


probe1=datatarget.probe;
probe2=datamovable.probe;

if(isa(probe1,'nirs.core.Probe1020'))
    probe1b=nirs.core.Probe;
    probe1b.optodes=probe1.optodes;
    probe1b.link=probe1.link;
    probe1=probe1b;
end
if(isa(probe2,'nirs.core.Probe1020'))
    probe2b=nirs.core.Probe;
    probe2b.optodes=probe2.optodes;
    probe2b.link=probe2.link;
    probe2=probe2b;
end



minX = min(min(probe1.optodes.X),min(probe2.optodes.X));
maxX = max(max(probe1.optodes.X),max(probe2.optodes.X));
dX = (maxX-minX)/3;
minY = min(min(probe1.optodes.Y),min(probe2.optodes.Y));
maxY = max(max(probe1.optodes.Y),max(probe2.optodes.Y));
dY = (maxY-minY)/3;

[X,Y,Z]=meshgrid([minX-dX:dX/30:maxX+dX],[minY-dY:dY/30:maxY+dY],[-10]);

mesh=nirs.core.Mesh;
mesh.nodes=[X(:) Y(:) Z(:)];

lambda=unique(probe1.link.type);


probe1.link=probe1.link(ismember(probe1.link.type,probe1.link.type(1)),:);
probe2.link=probe2.link(ismember(probe2.link.type,probe2.link.type(1)),:);

if(~isnumeric(lambda))
    probe1.link.type=ones(height(probe1.link),1)*808;
    probe2.link.type=ones(height(probe2.link),1)*808;
    lambda=808;
end

T=eye(3);

FwdModel=nirs.forward.ApproxSlab;
FwdModel.mesh=mesh;
FwdModel.prop=nirs.media.tissues.brain(.7,50,lambda);
FwdModel.Fm=0;
FwdModel2=FwdModel;
FwdModel.probe=probe1;
FwdModel2.probe=probe2;
J1=FwdModel.jacobian;
J1=J1.mua;
J2=FwdModel2.jacobian;
J2=J2.mua;

%% TODO- replace this with a better inverse solution

[iJ1 a1] = nirs.math.inverse(J1,datatargT);

[iJ2 a2] = nirs.math.inverse(J2,datamovT);


[optimizer, metric] = imregconfig('multimodal');
optimizer.MaximumIterations = 50;

% Find geometric transformation that maps moving to fixed.
T=zeros(3,3);
for j=1:size(a1,2)
    fixed=reshape(a1(:,j),size(X));
    moving=reshape(a2(:,j),size(X));
    tform = imregtform(moving, fixed, type, optimizer, metric);
    T=T+tform.T;
end
T=T/size(a1,2);
tform.T=T;

% if the two probes are the same then this works... but use a more general
% form here.
%datamoved=nirs.classify.util.apply_warp(datamovable,tform);

mesh2=mesh;
for j=1:2
    moving=reshape(mesh.nodes(:,j),size(X));
    mesh2.nodes(:,j) = reshape(imwarp(moving,tform.invert,'OutputView',imref2d(size(fixed))),[],1);
end
FwdModel2.mesh=mesh2;
J2=FwdModel2.jacobian;
J2=J2.mua;
a2=J1*iJ2;


n=length(datamovable.beta)/size(a2,1);
a2=kron(a2,eye(n));

datamoved=datamovable;
datamoved=datamoved.sorted();
datamoved.beta=a2*datamoved.beta;
datamoved.covb=a2*datamoved.covb*a2';





