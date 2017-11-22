function [VitE,mesh] = findVitaminE_T1wMRI(T1wfile,numMarkers,showplot)

if(nargin<2 || isempty(numMarkers))
    numMarkers=4;
end


if(nargin<3)
    showplot=true;
end

vol=load_untouch_nii(T1wfile);
vol2=medfilt3(vol.img,[5 5 5]);

[faces2,nodes2]=isosurface(vol2,50);

vol2(:,:,1:end/2)=0;

[faces,nodes]=isosurface(vol2,50);
nodes=nodes(:,[2 1 3]);
nodes(:,1)=nodes(:,1)*vol.hdr.dime.pixdim(1);
nodes(:,2)=nodes(:,2)*vol.hdr.dime.pixdim(2);
nodes(:,3)=nodes(:,3)*vol.hdr.dime.pixdim(3);
nodes(:,4)=1;
Tf=[vol.hdr.hist.srow_x; vol.hdr.hist.srow_y; vol.hdr.hist.srow_z];
Tf(4,4)=1;

nodes=(Tf*nodes')';
nodes(:,4)=[];


nodes2=nodes2(:,[2 1 3]);
nodes2(:,1)=nodes2(:,1)*vol.hdr.dime.pixdim(1);
nodes2(:,2)=nodes2(:,2)*vol.hdr.dime.pixdim(2);
nodes2(:,3)=nodes2(:,3)*vol.hdr.dime.pixdim(3);
nodes2(:,4)=1;
nodes2=(Tf*nodes2')';
nodes2(:,4)=[];

m=mean(nodes2,1);
nodes2=nodes2-ones(size(nodes2,1),1)*m;
nodes=nodes-ones(size(nodes,1),1)*m;



connected=zeros(size(nodes,1),1);
idx=1;
while(1)
    lst2=find(connected==0);
    if(isempty(lst2))
        break
    end
    
    disp(['Iteration ' num2str(idx) ' (' num2str(length(lst2)/length(connected)*100) '% remaining)']);
    [~,lst]=min(nodes(lst2,3));
    
    connected(lst2(lst))=idx;
    
    n=0; n1=length(find(connected==idx));
    while(n1>n)
        f=find(any(ismember(faces,find(connected==idx)),2));
        connected(faces(f,1))=idx;
        connected(faces(f,2))=idx;
        connected(faces(f,3))=idx;
        n=n1;
        n1=length(find(connected==idx));
        
    end
    
    idx=idx+1;
end


for i=1:idx-1; 
    cnt(i)=length(find(connected==i)); 
end;
[~,i]=sort(cnt,'descend');


VitE=struct;
connected2=zeros(size(connected));
for j=2:(1+numMarkers)
    connected2(ismember(connected,i(j)))=j-1;
    VitE.Name{j-1,1}='VitE';
    VitE.X(j-1,1)=mean(nodes(ismember(connected,i(j)),1));
    VitE.Y(j-1,1)=mean(nodes(ismember(connected,i(j)),2));
    VitE.Z(j-1,1)=mean(nodes(ismember(connected,i(j)),3));
    VitE.Type{j-1,1}='FID';
    VitE.Units{j-1,1}='mm';
end
VitE=struct2table(VitE);

if(showplot)
    figure; hold on;
    lst2=find(ismember(connected2,[1:numMarkers]));
    lst=find(any(ismember(faces,lst2),2));
    
    nirs.util.plotmesh(nodes2,faces2);
    nirs.util.plotmesh(nodes,faces(lst,:),connected2);
    caxis([0 numMarkers]);
    cm=jet(5);
    cm(1,:)=[.8 .8 .8];
    colormap(cm)
end

if(nargout==2)
    
    lst2=find(ismember(connected2,[1:numMarkers]));
    lst=find(ismember(nodes2,nodes(lst2,:),'rows'));
    nodes2(lst,:)=NaN;
    lst=find(any(ismember(faces2,lst),2));
    
    mesh=nirs.core.Mesh(nodes2,faces2);
    mesh=mesh.addfiducials(VitE);
end
