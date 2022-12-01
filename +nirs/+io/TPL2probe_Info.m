function probe1020=TPL2probe_Info(tplfile,mntfile)

ml=dlmread(tplfile);

tmp=load(mntfile);
mnt=tmp.(strtok(mntfile,'.'));

[i,j]=find(ml~=0);
for id=1:length(i)
    str=num2str(ml(i(id),j(id)));
    if(length(str)==3)
        str=['0' str];
    end
    sIdx = str2num(str(1:2));
    dIdx = str2num(str(3:4));
    meas(id,1)=sIdx;
    meas(id,2)=dIdx;
end
meas=unique(meas,'rows');  % sort by sources then detectors

% initial guess
uS=unique(meas(:,1));
for i=1:length(uS)
    lst=find(meas(:,1)==uS(i));
    srcPos(i,1)=mean(mnt.x(lst));
    srcPos(i,2)=mean(mnt.y(lst));
    srcPos(i,3)=0;
    
    srcPos3Db(i,1)=mean(mnt.pos_3d(1,lst)');
    srcPos3Db(i,2)=mean(mnt.pos_3d(2,lst)');
    srcPos3Db(i,3)=mean(mnt.pos_3d(3,lst)');
end

uD=unique(meas(:,2));
for i=1:length(uS)
    lst=find(meas(:,2)==uD(i));
    detPos(i,1)=mean(mnt.x(lst));
    detPos(i,2)=mean(mnt.y(lst));
    detPos(i,3)=0;
    
    detPos3Db(i,1)=mean(mnt.pos_3d(1,lst)');
    detPos3Db(i,2)=mean(mnt.pos_3d(2,lst)');
    detPos3Db(i,3)=mean(mnt.pos_3d(3,lst)');
end

link=table([meas(:,1);meas(:,1)],[meas(:,2);meas(:,2)],[repmat({'type-1'},size(meas,1),1); repmat({'type-2'},size(meas,1),1)],'VariableNames',{'source','detector','type'});


probe=nirs.core.Probe(srcPos,detPos,link);
scale=30/mean(probe.distances);   % the scale on this is really messed up, so rescale to the 30mm based on the hdr file

probe=nirs.core.Probe(srcPos*scale,detPos*scale,link);

xyz=[mnt.x mnt.y zeros(size(mnt.x))]*scale;
for i=1:length(mnt.clab)
    Name{i}=mnt.clab{i};
    Type{i}='FID-anchor';
    Units{i}='mm';
end
fid=table(Name',xyz(:,1),xyz(:,2),xyz(:,3),Type',Units',...
    'VariableNames',{'Name','X','Y','Z','Type','Units'});
probe.optodes=[probe.optodes; fid];

probe1020=nirs.util.registerprobe1020(probe);

Colin=nirs.registration.Colin27.mesh;
Colin(1).fiducials.Draw=ismember(Colin(1).fiducials.Name,mnt.clab);
probe1020=probe1020.register_mesh2probe(Colin);

probe1020.defaultdrawfcn='3D mesh (superior)';


% this mostly worked, but let's fix it a bit

[a,b]=ismember(mnt.clab,Colin(1).fiducials.Name);
lst=b(find(a));
p3D=[Colin(1).fiducials.X(lst) Colin(1).fiducials.Y(lst) Colin(1).fiducials.Z(lst)];

uS=unique(meas(:,1));
for i=1:length(uS)
    lst=find(meas(:,1)==uS(i));
    srcPos3D(i,:)=mean(p3D(lst,:),1);
end

uD=unique(meas(:,2));
for i=1:length(uD)
    lst=find(meas(:,2)==uD(i));
    detPos3D(i,:)=mean(p3D(lst,:),1);
end


lst=ismember(probe1020.optodes_registered.Type,'Source');
probe1020.optodes_registered.X(lst)=srcPos3D(:,1);
probe1020.optodes_registered.Y(lst)=srcPos3D(:,2);
probe1020.optodes_registered.Z(lst)=srcPos3D(:,3);

lst=ismember(probe1020.optodes_registered.Type,'Detector');
probe1020.optodes_registered.X(lst)=detPos3D(:,1);
probe1020.optodes_registered.Y(lst)=detPos3D(:,2);
probe1020.optodes_registered.Z(lst)=detPos3D(:,3);
    
probe1020 = nirs.registration.realign_symetric(probe1020);


return








