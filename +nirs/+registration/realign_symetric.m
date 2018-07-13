function probe1020=realign_symetric(probe1020)
% This function realigns a 10-20 probe to be symetric 

optodes=probe1020.optodes_registered;
optodes2D=probe1020.optodes;

pos3D=[optodes.X optodes.Y optodes.Z];
pos2D=[optodes2D.X optodes2D.Y optodes2D.Z];

pos2D = pos2D-ones(size(pos2D,1),1)*mean(pos2D,1);

% Find symetric pairs
d=squareform(pdist(abs(pos2D)));
[i,j]=find(d<5);
lst=find(i~=j);
pairs=[i(lst) j(lst)];

for idx=1:size(pairs,1)
    A=pos3D(pairs(idx,1),:);
    B=pos3D(pairs(idx,2),:);
    C=.5*(abs(A)+abs(B));
    A=C.*sign(A);
    B=C.*sign(B);
    pos3D(pairs(idx,1),:)=A;
    pos3D(pairs(idx,2),:)=B;   
end

% mesh=probe1020.getmesh;
% k=dsearchn(mesh(1).nodes,pos3D);
% pos3D=mesh(1).nodes(k,:);



probe1020.optodes_registered.X=pos3D(:,1);
probe1020.optodes_registered.Y=pos3D(:,2);
probe1020.optodes_registered.Z=pos3D(:,3);
