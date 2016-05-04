function probe1020=register2polhemus(probe,polfile)
% This function registers the probe to the polhemus text file

% Read the polhemus file
fid=fopen(polfile,'r');
C=textscan(fid,'%s %f %f %f');
fclose(fid);

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
        probe.optodes.X(id)=C{2}(lst(i));
        probe.optodes.Y(id)=C{3}(lst(i));
        probe.optodes.Z(id)=C{4}(lst(i));
    end
end
fid=table(Name',xyz(:,1),xyz(:,2),xyz(:,3),Type',Units',...
    'VariableNames',{'Name','X','Y','Z','Type','Units'});
% and concatinate it to the probe
probe.optodes=[probe.optodes; fid];

probe=nirs.util.probe_remove_unconnected(probe);
probe1020=nirs.util.registerprobe1020(probe);


return