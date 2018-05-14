function raw = register2BrainStorm(raw,folder)


% convert positions to scs space
% Based on: Nasion, left pre-auricular point (LPA), and right pre-auricular point (RPA).
%               Origin: Midway on the line joining LPA and RPA
%               Axis X: From the origin towards the Nasion (exactly through)
%               Axis Y: From the origin towards LPA in the plane defined by (NAS,LPA,RPA), and orthogonal to X axis
%               Axiz Z: From the origin towards the top of the head 


brainstorm nogui;
import_anatomy(0);

[sSubject, iSubject] = bst_get('Subject', 0);
head_mesh_fn = sSubject.Surface(sSubject.iScalp).FileName;
sHead = in_tess_bst(head_mesh_fn);

cort_mesh_fn = sSubject.Surface(sSubject.iCortex).FileName;
sCortex = in_tess_bst(cort_mesh_fn);

sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);

probe1020=raw.probe;
LE=probe1020.optodes_registered(ismember(probe1020.optodes_registered.Name,'LPA'),:);
LE=[LE.X LE.Y LE.Z];
RE=probe1020.optodes_registered(ismember(probe1020.optodes_registered.Name,'RPA'),:);
RE=[RE.X RE.Y RE.Z];
NA=probe1020.optodes_registered(ismember(probe1020.optodes_registered.Name,'NAS'),:);
NA=[NA.X NA.Y NA.Z];


mp=.5*( sMri.SCS.LPA+sMri.SCS.RPA);
a=sMri.SCS.NAS-.5*( sMri.SCS.LPA+sMri.SCS.RPA);  a=a/norm(a);
b=sMri.SCS.LPA-.5*( sMri.SCS.LPA+sMri.SCS.RPA); b=b/norm(b);
Targ = [sMri.SCS.NAS; sMri.SCS.LPA; sMri.SCS.RPA; cross(a,b)+mp];

Targ(:,1)=Targ(:,1)-mp(1);
Targ(:,2)=Targ(:,2)-mp(2);
Targ(:,3)=Targ(:,3)-mp(3);

a=NA-.5*(RE+LE);  a=a/norm(a);
b=LE-.5*(RE+LE);  b=b/norm(b);
Ref = [NA; LE; RE; cross(a,b)+.5*(RE+LE)];
Targ(:,4)=1;
Ref(:,4)=1;
R=Ref\Targ;
 
xyz=[probe1020.optodes_registered.X probe1020.optodes_registered.Y probe1020.optodes_registered.Z];
xyz(:,4)=1;
xyz=xyz*R;


probe1020.optodes_registered.X=xyz(:,2);
probe1020.optodes_registered.Y=xyz(:,1);
probe1020.optodes_registered.Z=xyz(:,3);

pts1020=nirs.util.list_1020pts('?');


LE=pts1020(ismember(pts1020.Name,'lpa'),:);
LE=[LE.X LE.Y LE.Z];
RE=pts1020(ismember(pts1020.Name,'rpa'),:);
RE=[RE.X RE.Y RE.Z];
NA=pts1020(ismember(pts1020.Name,'nas'),:);
NA=[NA.X NA.Y NA.Z];
CZ=pts1020(ismember(pts1020.Name,'Cz'),:);
CZ=[CZ.X CZ.Y CZ.Z];


LE2=probe1020.optodes_registered(ismember(lower(probe1020.optodes_registered.Name),'lpa'),:);
LE2=[LE2.X LE2.Y LE2.Z];
RE2=probe1020.optodes_registered(ismember(lower(probe1020.optodes_registered.Name),'rpa'),:);
RE2=[RE2.X RE2.Y RE2.Z];
NA2=probe1020.optodes_registered(ismember(lower(probe1020.optodes_registered.Name),'nas'),:);
NA2=[NA2.X NA2.Y NA2.Z];

a=NA2-.5*(RE2+LE2);  a=a/norm(a);
b=LE2-.5*(RE2+LE2);  b=b/norm(b);
v=[0:.1:200]'*[0 0 1]; 
v=v+ones(size(v,1),1)*(.5*(RE2+LE2));

[k,d]=dsearchn(v,sHead.Vertices*1000);
[~,i]=min(abs(d));
CZ2=v(k(i),:);

A=[LE; RE; NA; CZ];
B=[LE2; RE2; NA2; CZ2];
A(:,4)=1;
B(:,4)=1;
R=B\A;


v=sHead.Vertices*1000;
v(:,4)=1;
v=v*R;
xyz=[pts1020.X pts1020.Y pts1020.Z];
[TR,TT]=icp(v(:,1:3)',xyz(:,1:3)');
v=(inv(TR) * (v(:,1:3)' - TT*ones(1,size(v,1))))';
m(1)=nirs.core.Mesh(v,sHead.Faces);
m(1).transparency=.1;
m(1)=m(1).addfiducials(pts1020);    
m(1).fiducials.Draw(:)=false;

v=sCortex.Vertices*1000;
v(:,4)=1;
v=v*R;
v=(inv(TR) * (v(:,1:3)' - TT*ones(1,size(v,1))))';
m(2)=nirs.core.Mesh(v,sCortex.Faces);

xyz=[probe1020.optodes_registered.X probe1020.optodes_registered.Y probe1020.optodes_registered.Z];
xyz(:,4)=1;
xyz=xyz*R;
xyz=(inv(TR) * (xyz(:,1:3)' - TT*ones(1,size(xyz,1))))';

probe1020.optodes_registered.X=xyz(:,1);
probe1020.optodes_registered.Y=xyz(:,2);
probe1020.optodes_registered.Z=xyz(:,3);

probe1020=probe1020.register_mesh2probe(m,true);

probe1020.defaultdrawfcn='3D mesh (top)';

raw.probe=probe1020;



