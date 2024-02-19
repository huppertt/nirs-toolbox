function [roi,BrainMap] = convertNIFTI_to_roi(probe,BrainMap)

if(isstr(BrainMap))
    name=BrainMap;
    BrainMap=project2colin(BrainMap);
else
    name='custom';
end

mesh=probe.getmesh;
mesh=mesh(end);

fwdModel=nirs.forward.ApproxSlab;
fwdModel.probe=probe;
fwdModel.probe.link.type=[];
fwdModel.probe.link.type=repmat(808,height(fwdModel.probe.link),1);
fwdModel.mesh=mesh;
fwdModel.prop=nirs.media.tissues.brain(808);

J=fwdModel.jacobian;
weights=-J.mua*BrainMap;
%weights=weights/sum(weights);

roi=probe.link;
roi.type=[];
roi.weight=weights;
roi.Name=repmat({name},height(roi),1);

figure;
probe.draw3d;
hold on;
mesh.draw(BrainMap);


return

function BM = project2colin(niifile,layer)
% colin=nirs.registration.Colin27.mesh_V2;
% BM=project2colin('mindfulness.nii.gz');
% colin(end).draw(BM);

colin=nirs.registration.Colin27.mesh_V2;

if(nargin<2)
    layer=length(colin);
end

T =[ 0.9964    0.0178    0.0173   -0.0000
   -0.0169    0.9957   -0.0444   -0.0000
   -0.0151    0.0429    1.0215    0.0000
   -0.4232  -17.5022   11.6967    1.0000]';
n=load_nii(niifile);
R=[n.hdr.hist.srow_x; n.hdr.hist.srow_y; n.hdr.hist.srow_z; [0 0 0 1]];

n2=zeros(size(n.img));

Rot = inv(inv(T)*R);
nodes=colin(layer).nodes;
nodes(:,4)=1;
m3=nodes*Rot';
m3=max(m3,1);
for i=1:3
    m3(:,i)=min(m3(:,i),size(n2,i));
end

lst=sub2ind(size(n2),round(m3(:,1)),round(m3(:,2)),round(m3(:,3)));
n2(lst)=1;

BM=n.img(lst);
return