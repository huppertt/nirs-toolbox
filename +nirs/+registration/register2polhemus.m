function probe1020=register2polhemus(probe,polfile,rereg)
% This function registers the probe to the polhemus text file

% Read the polhemus file
fid=fopen(polfile,'r');
C=textscan(fid,'%s %f %f %f');
fclose(fid);

if(nargin<3)
    rereg=true;
end

for i=1:length(C{1})
    C{1}{i}(strfind(C{1}{i},':'))=[];
end
% Fix a few common notation issues
if(~isempty(find(ismember(lower(C{1}),'nz'))))
    C{1}{ismember(lower(C{1}),'nz')}='nas';
end

lab1020=nirs.util.list_1020pts;
lst=find(ismember(lower(C{1}),lower(lab1020.Name)));
lst2=find(ismember(lower(lab1020.Name),lower(C{1})));

% 
% xyz=[lab1020.X(lst2) lab1020.Y(lst2) lab1020.Z(lst2)];
% xyz2=[C{2}(lst) C{3}(lst) C{4}(lst)];
% [TR,TT]=icp(xyz',xyz2');



cnt=1;
for id=1:length(lst)
    Name{cnt}=C{1}{lst(id)};
    xyz(cnt,:)=[C{2}(lst(id)) C{3}(lst(id)) C{4}(lst(id))];
    Type{cnt}='FID-anchor';  % This is an anchor point
    Units{cnt}='mm';
    cnt=cnt+1;
end

probe=nirs.util.probe_remove_unconnected(probe);

fid=table(Name',xyz(:,1),xyz(:,2),xyz(:,3),Type',Units',...
    'VariableNames',{'Name','X','Y','Z','Type','Units'});
headsize=nirs.registration.getheadshape(fid);
probe1020 = nirs.core.Probe1020([],headsize);
probe1020.optodes=probe.optodes;
probe1020.optodes_registered=probe.optodes;
probe1020.optodes_registered(~ismember(probe1020.optodes_registered.Type,{'Source','Detector'}),:)=[];

probe1020.link=probe.link;

% Now, add the src-det points
lst=find(~ismember(lower(C{1}),lower(lab1020.Name)));
for i=1:length(lst)
    str=C{1}{lst(i)};
    if(strcmp(lower(str(1)),'s'))
        str2=['000' str(2:end)];
        str=['Source-' str2(end-3:end)];
    elseif(strcmp(lower(str(1)),'d'))
        str2=['000' str(2:end)];
        str=['Detector-' str2(end-3:end)];
    else
        str=[];
    end
    if(~isempty(str))
        id=find(ismember(probe.optodes.Name,str));
        probe1020.optodes_registered.X(id)=C{2}(lst(i));
        probe1020.optodes_registered.Y(id)=C{3}(lst(i));
        probe1020.optodes_registered.Z(id)=C{4}(lst(i));
        Name{cnt}=str;
        xyz(cnt,:)=[C{2}(lst(i)) C{3}(lst(i)) C{4}(lst(i))];
        Type{cnt}='FID-anchor';  % This is an anchor point
        Units{cnt}='mm';
        cnt=cnt+1;
        
    end
        
end

fid=table(Name',xyz(:,1),xyz(:,2),xyz(:,3),Type',Units',...
    'VariableNames',{'Name','X','Y','Z','Type','Units'});

% and concatinate it to the probe
probe1020.optodes_registered=[probe1020.optodes_registered; fid];
mesh=probe1020.getmesh;
T = nirs.registration.cp2tform(fid,mesh(1).fiducials);

mesh=nirs.registration.applytform(mesh,inv(T));
% mesh(1).fiducials.Draw(:)=false;
% fid.Draw(:)=true;
%mesh(1).fiducials=[mesh(1).fiducials; fid];

probe1020=probe1020.set_mesh(mesh,mesh(1).fiducials);


% 
% [xyz]=[probe1020.optodes_registered.X probe1020.optodes_registered.Y...
%     probe1020.optodes_registered.Z];
% xyz(:,4)=1;
% 
% if(isa(probe,'nirs.core.Probe1020'))
%     probe1020=nirs.registration.applytform(probe,T);
%     return
% end
% 
% xyz=xyz*T;
% probe1020.optodes_registered.X=xyz(:,1);
% probe1020.optodes_registered.Y=xyz(:,2);
% probe1020.optodes_registered.Z=xyz(:,3);





return