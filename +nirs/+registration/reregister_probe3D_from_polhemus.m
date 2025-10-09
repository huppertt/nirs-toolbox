function probeR=reregister_probe3D_from_polhemus(probe,poltextfile,rigid)

if(nargin<3)
    rigid=false;
end

optodes=probe.optodes_registered;

tblpol=readtable(poltextfile,'FileType','text');
tblpol.Properties.VariableNames={'Name','X','Y','Z'};
tblpol.Units=repmat({'mm'},height(tblpol),1);

tblpol.Name=strrep(tblpol.Name,'RPA','rpa');
tblpol.Name=strrep(tblpol.Name,'LPA','lpa');
tblpol.Name=strrep(tblpol.Name,'Nz','nas');

for i=1:height(tblpol)
    if(strcmp(lower(tblpol.Name{i}(1)),'s'))
        lst=max(find(double(tblpol.Name{i})>57))+1;
        idx=['0000' tblpol.Name{i}(lst:end)];
        tblpol.Name{i}=['Source-' idx(end-3:end)];
    elseif(strcmp(lower(tblpol.Name{i}(1)),'d'))
        lst=max(find(double(tblpol.Name{i})>57))+1;
        idx=['0000' tblpol.Name{i}(lst:end)];
        tblpol.Name{i}=['Detector-' idx(end-3:end)];
    end
end


BEM=probe.getmesh;
BEM(1).fiducials.Draw(:)=false;
BEM(1).fiducials(ismember(BEM(1).fiducials.Name,tblpol.Name),:).Draw(:)=true;

tbl1020=nirs.util.list_1020pts;
[~,tbl1020r] = nirs.registration.cp2tform(tbl1020,BEM(1).fiducials,true);

% This code takes the common named points in tblpol and finds the
% registration to those in the BEM mesh
T=nirs.registration.cp2tform(tblpol,tbl1020r,rigid);

% Then apply the transform so now tblpol(r) is registered to the BEM mesh
tblpolr=nirs.registration.applytform(tblpol,T);

% The points are potentially floating above/below the surface of the head,
% so refine the positions
xyz=[tblpolr.X tblpolr.Y tblpolr.Z];
xyz=nirs.registration.projectsurface(xyz,BEM(1).nodes);
tblpolr.X=xyz(:,1);
tblpolr.Y=xyz(:,2);
tblpolr.Z=xyz(:,3);

% Now that we have the digitizer point regsitered to the head, use those to
% bring all the other optodes points to the same space
[T,optodesnew]=nirs.registration.cp2tform(optodes,tblpolr,rigid);


[a,b]=ismember(optodesnew.Name,tblpolr.Name);

optodesnew(a,:).X=tblpolr(b(find(a)),:).X;
optodesnew(a,:).Y=tblpolr(b(find(a)),:).Y;
optodesnew(a,:).Z=tblpolr(b(find(a)),:).Z;

probeR=probe;
probeR.optodes_registered=optodesnew;
probeR2=probeR.register_mesh2probe(BEM);



% 
% %% Uncomment to draw alignment
% figure;
% BEM(1).draw; 
% hold on;
% xyz=[tblpolr.X tblpolr.Y tblpolr.Z];
% s=scatter3(xyz(:,1),xyz(:,2),xyz(:,3),'r','filled');
% xyz2=[optodesnew.X optodesnew.Y optodesnew.Z];
% s2=scatter3(xyz2(:,1),xyz2(:,2),xyz2(:,3),'g','filled');
% xyz3=[optodes.X optodes.Y optodes.Z];
% s3=scatter3(xyz3(:,1),xyz3(:,2),xyz3(:,3),'m','filled');
% legend({'','Atlas head points','Polhemus','Final Register','Before Polhemus'})