function probe2ProbeInfoMat(probe1020,filename)
% this saves a registered probe into the NIRx Probe_probeInfo.mat format

if(nargin<2)
    filename='probe_probeInfo.mat';
end

if(isempty(strfind(filename,'_probeInfo.mat')))
    filename=[strtok(filename,'.') '=_probeInfo.mat'];
end

lst=find(ismember(probe1020.optodes.Type,'Source'));
probeInfo.probes.coords_s2=[probe1020.optodes.X(lst) probe1020.optodes.Y(lst)];
probeInfo.probes.coords_s3=[probe1020.optodes_registered.X(lst) probe1020.optodes_registered.Y(lst) probe1020.optodes_registered.Z(lst)];

lst=find(ismember(probe1020.optodes.Type,'Detector'));
probeInfo.probes.coords_d2=[probe1020.optodes.X(lst) probe1020.optodes.Y(lst)];
probeInfo.probes.coords_d3=[probe1020.optodes_registered.X(lst) probe1020.optodes_registered.Y(lst) probe1020.optodes_registered.Z(lst)];

mesh=probe1020.getmesh;

probeInfo.geom.NIRxHead.ext1020sys.labels=mesh(1).fiducials.Name';
probeInfo.geom.NIRxHead.ext1020sys.coords3d=[mesh(1).fiducials.X mesh(1).fiducials.Y mesh(1).fiducials.Z];

for i=1:size(probeInfo.probes.coords_s3,1)
    probeInfo.geom.NIRxHead.ext1020sys.labels{end+1}=['Src-' num2str(i)];
    probeInfo.geom.NIRxHead.ext1020sys.coords3d(end+1,:)=probeInfo.probes.coords_s3(i,:);
end

for i=1:size(probeInfo.probes.coords_d3,1)
    probeInfo.geom.NIRxHead.ext1020sys.labels{end+1}=['Det-' num2str(i)];
    probeInfo.geom.NIRxHead.ext1020sys.coords3d(end+1,:)=probeInfo.probes.coords_d3(i,:);
end

probeInfo.probes.index_s(:,1)=1:size(probeInfo.probes.coords_s3,1);
probeInfo.probes.index_d(:,1)=1:size(probeInfo.probes.coords_d3,1);

[~,probeInfo.probes.index_s(:,2)]=ismember(probeInfo.probes.coords_s3,probeInfo.geom.NIRxHead.ext1020sys.coords3d,'rows');
[~,probeInfo.probes.index_d(:,2)]=ismember(probeInfo.probes.coords_d3,probeInfo.geom.NIRxHead.ext1020sys.coords3d,'rows');

probeInfo.geom.NIRxHead.mesh.nodes=mesh(1).nodes;
probeInfo.geom.NIRxHead.mesh.belems=mesh(1).faces;
probeInfo.geom.NIRxHead.mesh1.nodes=mesh(2).nodes;
probeInfo.geom.NIRxHead.mesh1.belems=mesh(2).faces;
probeInfo.geom.NIRxHead.mesh2.nodes=mesh(end).nodes;
probeInfo.geom.NIRxHead.mesh2.belems=mesh(end).faces;

probeInfo.newversion=true;

save(filename,'probeInfo');

