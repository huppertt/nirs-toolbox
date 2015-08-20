function probeOut = registerProbe2Mesh(mesh,probe)
% This function registers a probe with the fiducial points stored in a mesh

% First make sure we have at least 3 common labels to use for registration
NamesInMesh=mesh.fiducials.Name;
NamesInProbe=probe.optodes.Name;

CommonNames=intersect(lower(NamesInMesh),lower(NamesInProbe));
if(length(CommonNames)<3)
    error('At least 3 points needed for registration');
end
CommonNames=NamesInProbe(find(ismember(lower(NamesInProbe),CommonNames)));

%First, let's put the probe onto a sphere radius=head radius
nodes=mesh.nodes;
center=mean(nodes,1);
nodes=nodes-ones(size(nodes,1),1)*center;
HeadRadius = median(sqrt(sum(nodes.^2,2)));

ProbePos=[probe.optodes.X probe.optodes.Y probe.optodes.Z];
lstCM=find(ismember(probe.optodes.Units,{'cm'}));
ProbePos(lstCM,:)=ProbePos(lstCM,:)*10;

%Place the probe on a sphere
ProbePosSphere=[ProbePos(:,[1 2]) sign((HeadRadius.^2-sum(ProbePos(:,[1 2]).^2,2))).*real(sqrt(abs(HeadRadius.^2-sum(ProbePos(:,[1 2]).^2,2))))];


lstCommonInProbe=find(ismember(lower(NamesInProbe),lower(CommonNames)));
for idx=1:length(lstCommonInProbe)
    lstCommonInMesh(idx)=find(ismember(lower(NamesInMesh),lower(NamesInProbe(lstCommonInProbe(idx)))));
end

AnchorPos=[mesh.fiducials.X(lstCommonInMesh) mesh.fiducials.Y(lstCommonInMesh) mesh.fiducials.Z(lstCommonInMesh)];
AnchorPos=AnchorPos-ones(size(AnchorPos,1),1)*center;

%Deal with the anchors verses the attractors 
AnchorVectors=find(ismember(probe.optodes.Type(lstCommonInProbe),'FID-attractor'));
Anchors=find(ismember(probe.optodes.Type(lstCommonInProbe),'FID-anchor'));

dispflag=true;
if(dispflag)
    figure;
    mesh.transparency=.2;
    mesh.draw;
    hold on;
    s=scatter3(ProbePosSphere(:,1),ProbePosSphere(:,2),ProbePosSphere(:,3),'filled','r');
end



errPrev=inf;
for iter=1:10
    ProbeAnchors=ProbePosSphere(lstCommonInProbe,:);
    P1=ProbeAnchors(Anchors,:);
    P2=AnchorPos(Anchors,:);
    
    for i=1:length(Anchors)
        p1pos=P1(i,:);
        p2pos=P2(i,:);
        for j=1:length(AnchorVectors)
            p1v=ProbeAnchors(AnchorVectors(j),:)-p1pos;
            p2v=AnchorPos(AnchorVectors(j),:)-p2pos;
            p1v=p1v/norm(p1v);
            p2v=p2v/norm(p2v);
            P1=[P1; -p1v];
            P2=[P2; -p2v];
            
            
        end
    end
    P2(any(isnan(P1),2),:)=[];
    P1(any(isnan(P1),2),:)=[];
    
    Tform = P1\P2;
    %Tform=Tform/norm(Tform);
    ProbePosSphereReg=ProbePosSphere*Tform;
    [TR, TT] = icp(nodes',ProbePosSphereReg');
    ProbePosSphereReg=(TR*ProbePosSphereReg'+TT*ones(1,size(ProbePosSphereReg,1)))';
        
    if(dispflag)
        set(s,'XData',ProbePosSphereReg(:,1),'Ydata',ProbePosSphereReg(:,2),'zdata',ProbePosSphereReg(:,3));
    end
    ProbePosSphere=ProbePosSphereReg;
    
end

    for idx=1:size(ProbePosSphereReg,1)
        i=dsearchn(mesh.nodes, ProbePosSphereReg(idx,:));
        ProbePosSphereReg(idx,:)=mesh.nodes(i,:)-center;
    end


ProbePosSphereReg=ProbePosSphereReg+ones(size(ProbePosSphereReg,1),1)*center;

probeOut=probe;
probeOut.optodes.X=ProbePosSphereReg(:,1);
probeOut.optodes.Y=ProbePosSphereReg(:,2);
probeOut.optodes.Z=ProbePosSphereReg(:,3);

close;

return