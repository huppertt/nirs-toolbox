function probe1020=probe2T1wFiducials(T1wfile,digpts)

[v,mesh]=nirs.registration.findVitaminE_T1wMRI(T1wfile);

% first, lets sort the Vitamin E
lpa = digpts(ismember(digpts.Name,'lpa'),:);
rpa = digpts(ismember(digpts.Name,'rpa'),:);
nas = digpts(ismember(digpts.Name,'nas'),:);
vE1 = digpts(ismember(digpts.Name,'VitE'),:);

xyz1=[lpa.X lpa.Y lpa.Z; rpa.X rpa.Y rpa.Z;...
      nas.X nas.Y nas.Z; vE1.X vE1.Y vE1.Z];  


vE2 = mesh.fiducials(ismember(mesh.fiducials.Name,'VitE'),:);
xyz2=[vE2.X vE2.Y vE2.Z];

rr=xyz1(1:3,1:2)\[-1 0; 1 0; 0 1];
xyzr=xyz1(:,1:2)*rr;

th1=cart2pol(xyzr(:,1),xyzr(:,2));

th2=cart2pol(xyz2(:,1),xyz2(:,2));


[~,id1]=sort(th1(4:end));
xyz1(4:end,:)=xyz1(3+id1,:);

[~,id2]=sort(th2);
xyz2=xyz2(id2,:);


a=xyz1(4:end,:);
b=xyz2;

tform = pcregrigid(pointCloud(a),pointCloud(b));
a2=[digpts.X digpts.Y digpts.Z];
a2(:,4)=1;
%b = a*tform.T;
ar = a2*tform.T;
ar(:,4)=[];


% 
% [TR,TT]=icp(b',a');
% 
% a2=[digpts.X digpts.Y digpts.Z];
% ar=(TR*a2'+TT*ones(1,size(a2,1)))';
% 
[TR,TT]=icp(mesh.nodes',ar');
ar=(TR*ar'+TT*ones(1,size(ar,1)))';


lst=find(ismember(mesh.fiducials.Name,'VitE'));
for i=1:length(lst)
    mesh.fiducials.Name{lst(i)}=['VitE-' num2str(id2(i))];
end

k=dsearchn(mesh.nodes,ar);
ar=mesh.nodes(k,:);

digpts.X=ar(:,1);
digpts.Y=ar(:,2);
digpts.Z=ar(:,3);


lst=find(ismember(digpts.Type,{'Source','Detector'}));

probe1020 = nirs.core.Probe1020;
probe1020=probe1020.regsister_mesh2probe(mesh,true);
probe1020.optodes_registered=digpts(lst,:); 
 
