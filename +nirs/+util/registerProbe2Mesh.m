function probeOut = registerProbe2Mesh(mesh,probe,dispflag)

if(nargin<3)
    dispflag=false;
end

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
center=nanmean(nodes,1);
nodes=nodes-ones(size(nodes,1),1)*center;
HeadRadius = nanmedian(sqrt(sum(nodes.^2,2)));

ProbePos=[probe.optodes.X probe.optodes.Y probe.optodes.Z];
lstCM=find(ismember(probe.optodes.Units,{'cm'}));
ProbePos(lstCM,:)=ProbePos(lstCM,:)*10;

% Find symetric pairs;
sym=zeros(size(ProbePos,1));
for i=1:size(ProbePos,1)
    for j=i:size(ProbePos,1)
        if(norm(abs(ProbePos(i,:))-abs(ProbePos(j,:)))<5)
        sym(i,j)=1;
        sym(j,i)=1;
        end
    end
end

lst=find(ismember(probe.optodes.Type,{'Source','Detector','FID-anchor'}));

%Place the probe on a sphere
x=ProbePos(:,1); x=x-nanmean(x(lst));
y=ProbePos(:,2); y=y-nanmean(y(lst));
theta = atan(x/HeadRadius);
phi = atan(y/HeadRadius);
[ProbePosSphere(:,1),ProbePosSphere(:,2),ProbePosSphere(:,3)]=sph2cart(theta,phi,HeadRadius*ones(size(x)));

lstCommonInProbe=find(ismember(lower(NamesInProbe),lower(CommonNames)));
for idx=1:length(lstCommonInProbe)
    lstCommonInMesh(idx)=find(ismember(lower(NamesInMesh),lower(NamesInProbe(lstCommonInProbe(idx)))));
end

AnchorPos=[mesh.fiducials.X(lstCommonInMesh) mesh.fiducials.Y(lstCommonInMesh) mesh.fiducials.Z(lstCommonInMesh)];

%Deal with the anchors verses the attractors 
AnchorVectors=find(ismember(probe.optodes.Type(lstCommonInProbe),'FID-attractor'));
Anchors=find(ismember(probe.optodes.Type(lstCommonInProbe),'FID-anchor'));

% Initial registration
T = ProbePosSphere(lstCommonInProbe,:)\AnchorPos;
ProbePosSphere=ProbePosSphere*T;

ProbePosSphere=makesymetric(ProbePosSphere,sym);
ProbePosSphere = projectsurface(ProbePosSphere(:,1:3),mesh.nodes);

if(dispflag)
    figure;
    mesh.transparency=.2;
    mesh.draw;
    hold on;
    s=scatter3(ProbePosSphere(:,1),ProbePosSphere(:,2),ProbePosSphere(:,3),'filled','r');    
end


errPrev=inf;
for iter=1:1
    ProbeAnchors=ProbePosSphere(lstCommonInProbe,:);
 
    P1=ProbeAnchors(Anchors,:);
    P2=AnchorPos(Anchors,:);
    p1=ProbeAnchors(AnchorVectors,:);
    p2=AnchorPos(AnchorVectors,:);
    
    for j=1:length(AnchorVectors)
        k=dsearchn(P2(1:length(Anchors),:),p2(j,:));
        p1v = P1(k,:)-p1(j,:);
        p2v = P2(k,:)-p2(j,:);
        if(norm(p1v)~=0 & norm(p2v)~=0)
            p1v=p1v/norm(p1v);
            p2v=p2v/norm(p2v);
            P1=[P1; p1v-P1(k,:)];
            P2=[P2; p2v-P2(k,:)];
        end
    end
    Tform =P1\P2;
     
    ProbePosSphere=ProbePosSphere*Tform;
    ProbePosSphere=makesymetric(ProbePosSphere,sym);

    ProbePosSphere = projectsurface(ProbePosSphere,mesh.nodes);
    ProbePosSphere = pushdistances(ProbePosSphere,squareform(pdist(ProbePos)));
    ProbePosSphere=makesymetric(ProbePosSphere,sym);

    ProbePosSphere = projectsurface(ProbePosSphere,mesh.nodes);
    
    if(dispflag)
        s2=scatter3(AnchorPos(:,1),AnchorPos(:,2),AnchorPos(:,3),'filled','b');
        set(s,'XData',ProbePosSphere(:,1),'Ydata',ProbePosSphere(:,2),'zdata',ProbePosSphere(:,3));
        pause(2)
    end
     
end

P1=ProbeAnchors(Anchors,:);
P2=AnchorPos(Anchors,:);
p1=ProbeAnchors(AnchorVectors,:);
p2=AnchorPos(AnchorVectors,:);

for j=1:length(AnchorVectors)
    k=dsearchn(P2(1:length(Anchors),:),p2(j,:));
    p1v = P1(k,:)-p1(j,:);
    p2v = P2(k,:)-p2(j,:);
    if(norm(p1v)~=0 & norm(p2v)~=0)
        p1v=p1v/norm(p1v);
        p2v=p2v/norm(p2v);
        P1=[P1; p1v-P1(k,:)];
        P2=[P2; p2v-P2(k,:)];
    end
end
Tform =P1\P2;

ProbePosSphere=ProbePosSphere*Tform;
ProbePosSphere=makesymetric(ProbePosSphere,sym);
 
ProbePosSphere = projectsurface(ProbePosSphere,mesh.nodes);    
    
ProbePosSphere=makesymetric(ProbePosSphere,sym);
 
ProbePosSphere = projectsurface(ProbePosSphere,mesh.nodes);  

probeOut=probe;
probeOut.optodes.X=ProbePosSphere(:,1);
probeOut.optodes.Y=ProbePosSphere(:,2);
probeOut.optodes.Z=ProbePosSphere(:,3);

close;

return


function pos = pushdistances(pos,idealdist)
opt = optimoptions('lsqnonlin', 'MaxFunEvals', 1000,'Display','off');
 
mask=1*(idealdist(:)<45);

dx=zeros(length(pos(:)),1);
cost=@(dx)mask.*reshape(abs(squareform(pdist(pos+reshape(dx,size(pos))))-idealdist),[],1);
x=lsqnonlin(cost,dx,[],[],opt);
pos=pos+reshape(x,size(pos));

return


function pos = projectsurface(pos,surf)

com = nanmean(surf,1);

surf(any(isnan(surf),2),:)=[];
if(size(surf,1)>5000)
    surf=surf(randperm(size(surf,1),5000),:);
end
T=delaunayn(surf);

for idx=1:size(pos,1)
    disp(idx)
    vec = pos(idx,:)-com;
     c = [0:.1:2*norm(vec)];
    vec=vec/norm(vec);
    p=c'*vec+ones(length(c),1)*com;

   [k]=dsearchn(surf,T,p);
    d=sqrt(sum((surf(k,:)-p).^2,2));
    
    [~,i]=min(d);
    pos(idx,:)=p(i,:);
end


return


function ProbePosSphere=makesymetric(ProbePosSphere,sym);

for i=1:length(sym)
    lst=find(sym(i,:));
    p=nanmean(abs(ProbePosSphere(lst,:)),1);
    for j=1:length(lst)
        ProbePosSphere(lst(j),:)=p.*sign(ProbePosSphere(lst(j),:));
    end

end

return