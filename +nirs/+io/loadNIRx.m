function raw = loadNIRx(folder,registerprobe)
% This function loads NIRx data

% <subjid>.wl1  - wavelength #1 data
% <subjid>.wl2  - wavelength #2 data
% <subjid>_config.txt	- config file
% <subjid>.evt	- stimulus events	  (data taken from config file)
% <subjid>_probeInfo.mat - probe file
% <subjid>.tpl -topology file  (data taken from config file)

if(~isdir(folder))
    folder=fileparts(folder);
end

if(nargin<2)
    registerprobe=true;  % flag to also do 3D registration
end


file = dir(fullfile(folder,'*.hdr'));
if(isempty(file))
    file = dir(fullfile(folder,'*_config.json'));  % new format
    if(isempty(file))
        raw=[];
        return;
    else
        info = parsehdrJSON(fullfile(folder,file(1).name));
    end
else
    info = parsehdr(fullfile(folder,file(1).name));
end

%now read the probeInfo file and convert to a probe class
file = dir(fullfile(folder,'*robeInfo.mat'));
if(isempty(file)); raw=[]; return; end;

load(fullfile(folder,file(1).name)); % --> probeInfo

SrcPos = probeInfo.probes.coords_s2;
SrcPos(:,3)=0;

DetPos = probeInfo.probes.coords_d2;
DetPos(:,3)=0;


if(isfield(info,'ShortDetectors') && ~isempty(info.ShortDetectors) && ...
        info.ShortDetectors)
    info.S_D_Mask(:,end-info.ShortDetectors+1:end)=[];
    DetPos(info.Detectors-info.ShortDetectors+1:end,:)=[];
    useshortdistances=true;
else
    useshortdistances=false;
end
%
% info.S_D_Mask =[
%     1     1     0     1
%     1     1     1     0
%     1     0     0     0
%     1     0     0     0
%     0     0     1     0
%     0     0     1     0
%     0     0     0     1
%     0     0     0     1
%     ];


%this is a hack that allowed me to fix a study where the NIRx aquistion had
%been setup wrong and was missing some channels.  NIRx saves all data in
%the wl files (even those not part of the probe).  Creating a SDmask.mat
%file overwrites the info in the NIRx file
if(~isempty(dir(fullfile(folder,'*SDmask.mat'))))
    fi=dir(fullfile(folder,'*SDmask.mat'));
    load(fullfile(folder,fi(1).name));
    disp('loading mask from file');
    info.S_D_Mask=SD_mask;
end

[s,d]=find(info.S_D_Mask);


%% Fix added based on bug from Guilherme A. Zimeo Morais created an issue 2016-08-24
[s,a]=sort(s); d=d(a);

link=table(repmat(s,length(info.Wavelengths),1),...
    repmat(d,length(info.Wavelengths),1),...
    reshape(repmat(info.Wavelengths(:)',length(s),1),[],1),...
    'VariableNames',{'source','detector','type'});

probe = nirs.core.Probe(SrcPos,DetPos,link);


% if(isfield(info,'ShortDetIndex') && ~isempty(info.ShortDetIndex))
%      info.ShortDetectors=true;
% end



if(useshortdistances && info.ShortDetectors>0)
    
    d=(info.Detectors-info.ShortDetectors)+[1:info.ShortDetectors]';
    if(size(probeInfo.probes.coords_d2,1)==info.Detectors)
        s=dsearchn(probeInfo.probes.coords_s2,probeInfo.probes.coords_d2(d,:));
    else
        s=[1:info.ShortDetectors]';
    end
    ShortDetPos=SrcPos(s,:);
    %ShortSrcPos=SrcPos(s,:);
    
    Shortlink=table(repmat(s,length(info.Wavelengths),1),...
        repmat(d,length(info.Wavelengths),1),...
        reshape(repmat(info.Wavelengths(:)',length(s),1),[],1),...
        'VariableNames',{'source','detector','type'});
    probe.link=[probe.link table(repmat(false,height(probe.link),1),'VariableNames',{'ShortSeperation'})];
    Shortlink=[Shortlink table(repmat(true,height(Shortlink),1),'VariableNames',{'ShortSeperation'})];
    
    probe = nirs.core.Probe(SrcPos,[DetPos;ShortDetPos],[probe.link; Shortlink]);
    
end

lst2=[];
for j=1:length(probe.types)
    
    kk=find(ismember(probe.link.type,probe.types(j)));
    lst=find(ismember(info.SDkey(:,2:3),[probe.link.source(kk) probe.link.detector(kk)],'rows'));
    [~,i]=ismember(info.SDkey(lst,2:3),[probe.link.source(kk) probe.link.detector(kk)],'rows');
    lst2=[lst2; kk(i)];
end
probe.link=probe.link(lst2,:);


if(isfield(info,'ChanDis'))
    if(length(info.ChanDis)>length(probe.distances))
        info.ChanDis=info.ChanDis(1:length(probe.distances));
    end
    if(length(info.ChanDis)<length(probe.distances))
        info.ChanDis=reshape(repmat(info.ChanDis,1,length(info.Wavelengths)),[],1);
    end
    if(length(probe.distances)<length(info.ChanDis))
    info.ChanDis=info.ChanDis(1:probe.distances);
    end
    
    % Not sure why the units on the 2D probe in the NIRx file are so off
    l=info.ChanDis(:)./probe.distances(1:length(info.ChanDis(:)));
    l(l==Inf)=[];  % remove the short distance ones
    scale=mean(l);
    probe.fixeddistances=info.ChanDis(:);
else
    scale=1;
end


probe.optodes.X=scale*probe.optodes.X;
probe.optodes.Y=scale*probe.optodes.Y;
probe.optodes.Z=scale*probe.optodes.Z;

probe.link=sortrows(probe.link,{'type','source','detector'});

%% Now, let's get the data
file = rdir(fullfile(folder,'*.nirs' ));
if(~isempty(file))
    raw=nirs.io.loadDotNirs(file(1).name,true);
    link=raw.probe.link;
    probe.link=link;
    probe.optodes=raw.probe.optodes;
    XYZprobe3D=zeros(0,3);
    if(exist(fullfile(fileparts(file(1).name),'digpts.txt')))
        fil=fullfile(fileparts(file(1).name),'digpts.txt');
        fid=fopen(fil,'r');
        Name={}; Units={}; Type={};
        while(1)
            line=fgetl(fid);
            if(line==-1)
                break
            end
            XYZprobe3D(end+1,:)=str2num(line(strfind(line,':')+1:end));
            Name{end+1,1}=line(1:strfind(line,':')-1);
            Units{end+1,1}='mm';
            Type{end+1,1}='FID-anchor';
            if(strcmp(Name{end}(1),'S'))
                str=['0000' Name{end}(2:end)];
                Name{end}=['Source-' str(end-3:end)];
                Type{end}='Source';
            elseif(strcmp(Name{end}(1),'D'))
                str=['0000' Name{end}(2:end)];
                Name{end}=['Detector-' str(end-3:end)];
                Type{end}='Detector';
            end
            if(strcmp(lower(Name{end}),'nz'))
                Name{end}='Nas';
            end
            
        end
        X=XYZprobe3D(:,1);
        Y=XYZprobe3D(:,2);
        Z=XYZprobe3D(:,3);
        probe.optodes=table(Name,X,Y,Z,Type,Units);
        
        hs=nirs.registration.getheadshape(probe.optodes);
        t1020=nirs.util.list_1020pts('?',hs);
        [T probe.optodes] =  nirs.registration.cp2tform(probe.optodes,t1020);
        
        
    else
        
        com= [0 0 0];  %from HOMER-AtlasViewer
        
        XYZprobe3D=[raw.probe.optodes.X raw.probe.optodes.Y raw.probe.optodes.Z];
        
        XYZprobe3D(:,2)=XYZprobe3D(:,2)-ones(size(XYZprobe3D,1),1)*(com(2)+45);
        
        XYZprobe3D(:,1)=XYZprobe3D(:,1)-ones(size(XYZprobe3D,1),1)*(com(1));
        XYZprobe3D(:,2)=XYZprobe3D(:,2)-ones(size(XYZprobe3D,1),1)*(com(2));
        XYZprobe3D(:,3)=XYZprobe3D(:,3)-ones(size(XYZprobe3D,1),1)*(com(3));
        
        XYZprobe3D(:,2)=-XYZprobe3D(:,2);
        XYZprobe3D(:,[2 3])=XYZprobe3D(:,[3 2]);
        
        XYZprobe3D(:,1)=-XYZprobe3D(:,1);
        raw.probe.optodes.X=XYZprobe3D(:,1);
        raw.probe.optodes.Y=XYZprobe3D(:,2);
        raw.probe.optodes.Z=XYZprobe3D(:,3);
        
    end
    
    newVer=true;
else
    newVer=false;
    % read the standard wl files
    raw = nirs.core.Data();
    
    
    kk=find(ismember(probe.link.type,probe.link.type(1)));
    lst=find(ismember(info.SDkey(:,2:3),[probe.link.source(kk) probe.link.detector(kk)],'rows'));
    for idx=1:length(info.Wavelengths)
        file = dir(fullfile(folder,['*.wl' num2str(idx)]));
        
        if(file(1).bytes==0)
            error('zero byte file');
        end
        
        d = dlmread(fullfile(folder,file(1).name));
        raw.data=[raw.data d(:,lst)];
    end
    
    raw.time=[0:size(raw.data,1)-1]/info.SamplingRate;
end

if(isfield(info,'FileName'))
    raw.description=info.FileName;
else
    raw.description=folder;
end

% Add the demographics info
file = dir(fullfile(folder,'*.inf'));
if(~isempty(file))
    demoinfo = parsehdr(fullfile(folder,file(1).name));
else
    file = dir(fullfile(folder,'*_description.json'));
    if(~isempty(file))
        demoinfo = nirs.io.loadjson(fullfile(folder,file(1).name));
    end
end

demo=Dictionary;
flds=fields(demoinfo);
for idx=1:length(flds)
    demo(flds{idx})=demoinfo.(flds{idx});
end
raw.demographics=demo;

% There is a slight error in the NIRx files for hyperscanning that I am
% going to exploit to add this info.  For hyperscanning files, the
% info.S_D_Mask covers only 1 subject but the info.SD_Key field is the full
% (both subjects) length.
if(isfield(info,'ChanDis'))
    if(length(info.ChanDis)==length(s)*2)
        raw.demographics('hyperscan')=info.FileName;
    end
end


% Now add stimulus info
if(isfield(info,'Events') && ~isempty(info.Events))
    stimName = unique(info.Events(:,2));
    stimulus=Dictionary();
    for idx=1:length(stimName)
        s = nirs.design.StimulusEvents();
        s.name=['channel_' num2str(stimName(idx))];
        s.onset=info.Events(find(info.Events(:,2)==stimName(idx)),1);
        s.dur=ones(size(s.onset));
        s.amp=ones(size(s.onset));
        stimulus(s.name)=s;
    end
    raw.stimulus=stimulus;
end


if(~newVer)
    
    if(isfield(probeInfo,'geom'))
        
        %% Get the 10-20 fiducial points
        % Create the head mesh from the file
        Name=probeInfo.geom.NIRxHead.ext1020sys.labels';
        X=probeInfo.geom.NIRxHead.ext1020sys.coords3d(:,1);
        Y=probeInfo.geom.NIRxHead.ext1020sys.coords3d(:,2);
        Z=probeInfo.geom.NIRxHead.ext1020sys.coords3d(:,3);
        Type=repmat('10-20',length(Name),1);
        Units=repmat('mm',length(Name),1);
        Draw=repmat(true,length(Name),1);
        fid_1020 = table(Name,X,Y,Z,Type,Units,Draw);
        
        % Add FID locations to the 2D probe
        for i=1:size(probeInfo.probes.index_s,1)
            if(probeInfo.probes.index_s(i,2)>0)
                name{i,1}=fid_1020.Name(probeInfo.probes.index_s(i,2));
            else
                name{i,1}=['MNI:[' num2str(probeInfo.probes.coords_s3(i,1)) ','...
                    num2str(probeInfo.probes.coords_s3(i,2)) ',',...
                    num2str(probeInfo.probes.coords_s3(i,3)) ']'];
            end
        end
        
        fidS=table(name,...
            probe.srcPos(probeInfo.probes.index_s(:,1),1),...
            probe.srcPos(probeInfo.probes.index_s(:,1),2),...
            probe.srcPos(probeInfo.probes.index_s(:,1),3),...
            repmat({'FID-anchor'},length(probeInfo.probes.index_s),1),repmat({'mm'},length(probeInfo.probes.index_s),1),...
            'VariableNames',{'Name','X','Y','Z','Type','Units'});
        
        for i=1:size(probeInfo.probes.index_d,1)
            if(probeInfo.probes.index_d(i,2)>0)
                name2{i,1}=fid_1020.Name(probeInfo.probes.index_d(i,2));
            else
                name2{i,1}=['MNI:[' num2str(probeInfo.probes.coords_d3(i,1)) ','...
                    num2str(probeInfo.probes.coords_d3(i,2)) ',',...
                    num2str(probeInfo.probes.coords_d3(i,3)) ']'];
            end
        end
        
        for i=1:length(name2)
            if(iscell(name2{i}))
                name2{i}=name2{i}{1};
            end
        end
        
        for i=1:length(name)
            if(iscell(name{i}))
                name{i}=name{i}{1};
            end
        end
        
        fidD=table(name2,...
            probe.detPos(probeInfo.probes.index_d(:,1),1),...
            probe.detPos(probeInfo.probes.index_d(:,1),2),...
            probe.detPos(probeInfo.probes.index_d(:,1),3),...
            repmat({'FID-anchor'},length(probeInfo.probes.index_d),1),repmat({'mm'},length(probeInfo.probes.index_d),1),...
            'VariableNames',{'Name','X','Y','Z','Type','Units'});
        
        XYZprobe3D = [probeInfo.probes.coords_s3; probeInfo.probes.coords_d3; probeInfo.probes.coords_s3; probeInfo.probes.coords_d3];
        
        % and concatinate it to the probe
        probe.optodes=[probe.optodes; fidS; fidD];
        probeInfo.newversion=false;
        probeInfo.probes.index_s(:,1)=[];
        probeInfo.probes.index_d(:,1)=[];
    else
        
        load NIRxGeomHead;
        probeInfo.geom=geom;
        probeInfo.newversion=true;
        Name=probeInfo.geom.NIRxHead.ext1020sys.labels';
        X=probeInfo.geom.NIRxHead.ext1020sys.coords3d(:,1);
        Y=probeInfo.geom.NIRxHead.ext1020sys.coords3d(:,2);
        Z=probeInfo.geom.NIRxHead.ext1020sys.coords3d(:,3);
        Type=repmat('10-20',length(Name),1);
        Units=repmat('mm',length(Name),1);
        Draw=repmat(true,length(Name),1);
        fid_1020 = table(Name,X,Y,Z,Type,Units,Draw);
        
        % Add FID locations to the 2D probe
        fidS=table(probeInfo.probes.labels_s',scale*probeInfo.probes.coords_s2(:,1),...
            scale* probeInfo.probes.coords_s2(:,2),...
            scale*probeInfo.probes.coords_s2(:,2)*0,...
            repmat({'FID-anchor'},length(probeInfo.probes.coords_s3),1),repmat({'mm'},length(probeInfo.probes.coords_s3),1),...
            'VariableNames',{'Name','X','Y','Z','Type','Units'});
        fidD=table(probeInfo.probes.labels_d',scale*probeInfo.probes.coords_d2(:,1),...
            scale*probeInfo.probes.coords_d2(:,2),...
            scale*probeInfo.probes.coords_d2(:,2)*0,...
            repmat({'FID-anchor'},length(probeInfo.probes.coords_d3),1),repmat({'mm'},length(probeInfo.probes.coords_d3),1),...
            'VariableNames',{'Name','X','Y','Z','Type','Units'});
        
        fidS2=table(probeInfo.probes.labels_s',scale*probeInfo.probes.coords_s2(:,1),...
            scale* probeInfo.probes.coords_s2(:,2),...
            scale*probeInfo.probes.coords_s2(:,2)*0,...
            repmat({'Source'},length(probeInfo.probes.coords_s3),1),repmat({'mm'},length(probeInfo.probes.coords_s3),1),...
            'VariableNames',{'Name','X','Y','Z','Type','Units'});
        fidD2=table(probeInfo.probes.labels_d',scale*probeInfo.probes.coords_d2(:,1),...
            scale*probeInfo.probes.coords_d2(:,2),...
            scale*probeInfo.probes.coords_d2(:,2)*0,...
            repmat({'Detector'},length(probeInfo.probes.coords_d3),1),repmat({'mm'},length(probeInfo.probes.coords_d3),1),...
            'VariableNames',{'Name','X','Y','Z','Type','Units'});
        
        
        if(~newVer)
            XYZprobe3D = [probeInfo.probes.coords_s3; probeInfo.probes.coords_d3; probeInfo.probes.coords_s3; probeInfo.probes.coords_d3];
        end
        
        
        
        % and concatinate it to the probe
        probe.optodes=[probe.optodes(:,1)  [fidS2(:,2:end); fidD2(:,2:end)]; fidS; fidD];
        probe.optodes(ismember(probe.optodes.Name,''),:)=[];
        [log,idx] = ismember(probeInfo.probes.labels_s,probeInfo.geom.NIRxHead.ext1020sys.labels);
        probeInfo.probes.index_s = idx(log);
        [log,idx] = ismember(probeInfo.probes.labels_d,probeInfo.geom.NIRxHead.ext1020sys.labels);
        probeInfo.probes.index_d = idx(log);
        
    end
    
    
    
    if(registerprobe)
        
        for i=1:height(probe.optodes)
            if(iscell(probe.optodes.Name{i}))
                probe.optodes.Name{i}=probe.optodes.Name{i}{1};
            end
        end
        
        
        lst=find(ismember(probe.optodes.Type,'FID-anchor') & ~ismember(probe.optodes.Name,fid_1020.Name));
        if(isempty(lst))
            probe1020=nirs.util.registerprobe1020(probe);
        else
            
            t1020=nirs.util.list_1020pts('?');
            T=nirs.registration.cp2tform(fid_1020,t1020);
            XYZprobe3D(:,4)=1;
            XYZprobe3D=XYZprobe3D*T;
            try;
                extrapts=probe.optodes(lst,:);
                extrapts.X=XYZprobe3D(lst,1); extrapts.Y=XYZprobe3D(lst,2); extrapts.Z=XYZprobe3D(lst,3);
                probe1020=nirs.util.registerprobe1020(probe,[],[],extrapts);
            catch
                probe1020=nirs.util.registerprobe1020(probe);
                lst=[];
            end
        end
        
        if(isfield(probeInfo,'geom'))
            % old NIRx data format
            if(~isempty(lst))
                probeInfo.geom.NIRxHead.ext1020sys.coords3d=[probeInfo.geom.NIRxHead.ext1020sys.coords3d; XYZprobe3D(lst,1:3)];
                probeInfo.geom.NIRxHead.ext1020sys.labels={probeInfo.geom.NIRxHead.ext1020sys.labels{:} extrapts.Name{:}};
            end
            BEM(1)=nirs.core.Mesh(probeInfo.geom.NIRxHead.mesh.nodes(:,end-2:end),...
                probeInfo.geom.NIRxHead.mesh.belems(:,end-2:end),[]);
            %BEM(1)=reducemesh(BEM(1),.25);
            BEM(1).transparency=.2;
            BEM(1).fiducials=fid_1020;
            
            BEM(2)=nirs.core.Mesh(probeInfo.geom.NIRxHead.mesh1.nodes(:,end-2:end),...
                probeInfo.geom.NIRxHead.mesh1.belems(:,end-2:end),[]);
            %BEM(2)=reducemesh(BEM(2),.25);
            BEM(2).transparency=.2;
            
            BEM(3)=nirs.core.Mesh(probeInfo.geom.NIRxHead.mesh2.nodes(:,end-2:end),...
                probeInfo.geom.NIRxHead.mesh2.belems(:,end-2:end),[]);
            %BEM(3)=reducemesh(BEM(3),.25);
            BEM(3).transparency=1;
        else
            %new NIRx data format
            C27=nirs.registration.Colin27.BEM;
            BEM=C27.mesh;
        end
        
        % This will allow NIRFAST to directly use the info for the BEM model
        lambda=unique(probe.link.type);
        prop{1} = nirs.media.tissues.skin(lambda);
        prop{2} = nirs.media.tissues.bone(lambda);
        prop{3} = nirs.media.tissues.water(lambda);
        prop{4} = nirs.media.tissues.brain(lambda,0.7, 50);
        
        fwdBEM=nirs.forward.NirfastBEM;
        fwdBEM.mesh=BEM;
        fwdBEM.prop  = prop;
        
        probe1020=probe1020.register_mesh2probe(fwdBEM.mesh);
        probe1020.opticalproperties=prop;
        
        m=probe1020.getmesh;
        fid_1020=m(1).fiducials;
        
        if(useshortdistances)
            
            try
                probeInfo.geom.NIRxHead.ext1020sys.labels{find(ismember(probeInfo.geom.NIRxHead.ext1020sys.labels,{'FAF1'}))}='AFF1';
                probeInfo.geom.NIRxHead.ext1020sys.labels{find(ismember(probeInfo.geom.NIRxHead.ext1020sys.labels,{'FAF2'}))}='AFF2';
                probeInfo.geom.NIRxHead.ext1020sys.labels{find(ismember(probeInfo.geom.NIRxHead.ext1020sys.labels,{'FAF5'}))}='AFF5';
                probeInfo.geom.NIRxHead.ext1020sys.labels{find(ismember(probeInfo.geom.NIRxHead.ext1020sys.labels,{'FAF6'}))}='AFF6';
            end
            
            if(size(probeInfo.probes.index_s,1)>size(probeInfo.probes.index_s,2))
                probeInfo.probes.index_s=probeInfo.probes.index_s';
            end
            
            if(size(probeInfo.probes.index_d,1)>size(probeInfo.probes.index_d,2))
                probeInfo.probes.index_d=probeInfo.probes.index_d';
            end
            
            if(~isfield(probeInfo.probes,'index_s') || isempty(probeInfo.probes.index_s))
                index_s= find(ismember(probeInfo.geom.NIRxHead.ext1020sys.labels,probeInfo.probes.labels_s))';
                index_d= find(ismember(probeInfo.geom.NIRxHead.ext1020sys.labels,probeInfo.probes.labels_d))';
            else
                index_s=probeInfo.probes.index_s;
                index_d=probeInfo.probes.index_d;
                if(length(index_d)~=info.Detectors)
                    index_d=[index_d(1:end-1) index_s];
                end
            end
            
            lst=[index_s index_d...
                index_s(1:info.ShortDetectors)];
            labels={probeInfo.geom.NIRxHead.ext1020sys.labels{lst}};
            [~,lst2]=ismember(labels,fid_1020.Name);
            lst=[index_s index_d];
            labels={probeInfo.geom.NIRxHead.ext1020sys.labels{lst}};
            [~,lst3]=ismember(labels,fid_1020.Name);
            lst2=[lst2 lst3];
            
        else
            lst=[probeInfo.probes.index_s(:); probeInfo.probes.index_d(:)];
            if(~any(lst==0))
                labels={probeInfo.geom.NIRxHead.ext1020sys.labels{lst}};
                [~,lst2]=ismember(labels,fid_1020.Name);
                if(length(lst)*2==info.Detectors+info.Sources)
                    lst2=[lst2 lst2];
                end
            else
                lst2=[];
            end
        end
        
        if(~probeInfo.newversion)
            if(~isempty(lst2))
                XYZ=[fid_1020.X(lst2) fid_1020.Y(lst2) fid_1020.Z(lst2)];
                
                %
                %     SrcPos3D = probeInfo.probes.coords_s3;
                %     DetPos3D = probeInfo.probes.coords_d3;
                %     FID3D = [probeInfo.geom.NIRxHead.ext1020sys.coords3d(probeInfo.probes.index_s(:,2),:);...
                %         probeInfo.geom.NIRxHead.ext1020sys.coords3d(probeInfo.probes.index_d(:,2),:)];
                %
                %     XYZ=[SrcPos3D; DetPos3D; FID3D];
                %
                %     fidPts=probe1020.optodes_registered(ismember(probe1020.optodes_registered.Type,'FID-anchor'),:);
                %     XYZ_reg=[fidPts.X fidPts.Y fidPts.Z];x
                %     XYZ_reg(:,4)=1;
                %     FID3D(:,4)=1;
                %     XYZ(:,4)=1;
                %     R=FID3D\XYZ_reg;
                %
                %     XYZ=XYZ*R;
                
                probe1020.optodes_registered=probe1020.optodes;
                probe1020.optodes_registered.X(1:length(lst2))=XYZ(:,1);
                probe1020.optodes_registered.Y(1:length(lst2))=XYZ(:,2);
                probe1020.optodes_registered.Z(1:length(lst2))=XYZ(:,3);
            else
                % keep as is.
                XYZ=XYZprobe3D;
                probe1020.optodes_registered=probe1020.optodes;
                probe1020.optodes_registered.X(1:length(XYZ))=XYZ(:,1);
                probe1020.optodes_registered.Y(1:length(XYZ))=XYZ(:,2);
                probe1020.optodes_registered.Z(1:length(XYZ))=XYZ(:,3);
            end
        end
        if(newVer)
            XYZ=XYZprobe3D;
            probe1020.optodes_registered=probe1020.optodes;
            probe1020.optodes_registered.X(1:length(XYZ))=XYZ(:,1);
            probe1020.optodes_registered.Y(1:length(XYZ))=XYZ(:,2);
            probe1020.optodes_registered.Z(1:length(XYZ))=XYZ(:,3);
        end
        
        ll=[];
        for i=1:height(probe1020.optodes_registered)
            
            if(isempty(probe1020.optodes_registered.Name{i}))
                ll=[ll i];
            end
        end
        probe1020.optodes_registered(ll,:)=[];
        
        probe1020.link=probe.link;
        probe1020.fixeddistances=probe.distances;
        raw.probe=probe1020;
    else
        raw.probe=probe;
    end
else
    probe1020=nirs.core.Probe1020;
    probe1020.link=probe.link;
    probe1020.optodes_registered=probe.optodes;
    probe1020.optodes=probe.optodes;
    
    C27=nirs.registration.Colin27.BEM;
    BEM=C27.mesh;
    
    lst=ismember(probe1020.optodes.Type,{'Source','Detector'});
    
    probe1020.optodes.Z(:)=0;
    probe1020.optodes.X(lst)=[probeInfo.probes.coords_s2(:,1); probeInfo.probes.coords_d2(:,1)];
    probe1020.optodes.Y(lst)=[probeInfo.probes.coords_s2(:,2); probeInfo.probes.coords_d2(:,2)];
    
    s=raw.probe.distances\probe1020.distances;
    probe1020.optodes.X=probe1020.optodes.X/s;
    probe1020.optodes.Y=probe1020.optodes.Y/s;
    probe1020.fixeddistances=raw.probe.distances;
    
    if(~newVer)
        
        X=10*[probeInfo.probes.coords_d3(:,1); probeInfo.probes.coords_s3(:,1)];
        Y=10*[probeInfo.probes.coords_d3(:,2); probeInfo.probes.coords_s3(:,2)];
        Z=10*[probeInfo.probes.coords_d3(:,3); probeInfo.probes.coords_s3(:,3)];
        
        if(~isfield(probeInfo.probes,'labels_d'))
            for i=1:probeInfo.probes.nDetector0
                probeInfo.probes.labels_d{i}=['Detector-' num2str(i)];
            end
            
        end
        if(~isfield(probeInfo.probes,'labels_s'))
            for i=1:probeInfo.probes.nSource0
                probeInfo.probes.labels_s{i}=['Source-' num2str(i)];
            end
            
        end
        Name={probeInfo.probes.labels_d{:} probeInfo.probes.labels_s{:}}';
        Units=repmat({'mm'},size(Name));
        fid_1020=table(Name,X,Y,Z,Units);
    else
        fid_1020=probe.optodes;
    end
    
    t1020=nirs.util.list_1020pts('?');
    [T,a]=nirs.registration.cp2tform(t1020,fid_1020);
    BEM=nirs.registration.rotatemesh(BEM,T);
    
    % This will allow NIRFAST to directly use the info for the BEM model
    lambda=unique(probe1020.link.type);
    prop{1} = nirs.media.tissues.skin(lambda);
    prop{2} = nirs.media.tissues.bone(lambda);
    prop{3} = nirs.media.tissues.water(lambda);
    prop{4} = nirs.media.tissues.brain(lambda,0.7, 50);
    
    fwdBEM=nirs.forward.NirfastBEM;
    fwdBEM.mesh=BEM;
    fwdBEM.prop  = prop;
    
    probe1020=probe1020.register_mesh2probe(fwdBEM.mesh,true);
    probe1020.opticalproperties=prop;
    
    raw.probe=probe1020;
    
end

if(exist(fullfile(folder,'stimulus.mat')))
    load(fullfile(folder,'stimulus.mat'))
    raw.auxillary('stim')=raw.stimulus;
    raw.stimulus=stimulus;
    disp('loading stim-events from file');
end
if(~isempty(rdir(fullfile(folder,['*_nirsInfo.mat']))))
    f=rdir(fullfile(folder,['*_nirsInfo.mat']));
    info=load(f(1).name);
    if(isfield(info.nirsInfo,'eventInfo'))
        for i=1:length(info.nirsInfo.eventInfo.events.names)
            st=nirs.design.StimulusEvents;
            st.name=info.nirsInfo.eventInfo.events.names{i};
            st.onset=info.nirsInfo.eventInfo.events.onsets{i};
            st.dur=info.nirsInfo.eventInfo.events.durations{i};
            st.amp=ones(size(st.onset));
            raw.stimulus(st.name)=st;
        end
    end
end
end


function info =parsehdrJSON(file)
% This sub-routine parses the NIRx header info

info=nirs.io.loadjson(file);

info.S_D_Mask=[];
for i=1:length(info.channel_mask)
    for j=1:length(info.channel_mask{i});
        info.S_D_Mask(i,j)=str2num(info.channel_mask{i}(j));
    end
end

info.det_Mask=[];
for i=1:length(info.det_split)
    for j=1:length(info.det_split{i});
        info.det_Mask(i,j)=str2num(info.det_split{i}(j));
    end
end
lst=find(sum(info.det_Mask)==8);
if(~isempty(lst))
    info.ShortDetIndex=find(info.det_Mask(:,lst));
    info.ShortDetectors=length(info.ShortDetIndex);
else
    info.ShortDetIndex=[];
end
info.Detectors=length(sum(info.det_Mask));
info.Wavelengths=[760 850];
info.Sources=length(info.drv_amplitudes);

info.SDkey=[1:info.Sources*info.Detectors;repmat(1:info.Sources,1,info.Detectors); repmat(1:info.Detectors,1,info.Sources)]';

end

function info =parsehdr(file)
% This sub-routine parses the NIRx header info

info=struct;
fid=fopen(file,'r');

while(1)
    line=fgetl(fid);
    if(~isstr(line))
        break
    end
    if(~isempty(strfind(line,'=')))
        if(isempty(strfind(line,'="#')))
            fld = line(1:strfind(line,'=')-1);
            val = line(strfind(line,'=')+1:end);
            if(isempty(strfind(val,'"')))
                val=str2num(val);
            else
                val(strfind(val,'"'))=[];
            end
        else
            % matrix- read until next #
            fld = line(1:strfind(line,'=')-1);
            cnt=1; val=[];
            while(1)
                line=fgetl(fid);
                if(isempty(strfind(line(1),'#')))
                    val(cnt,:)=str2num(line);
                    cnt=cnt+1;
                else
                    break;
                end
            end
        end
        fld(strfind(fld,'-'))='_';
        fld(strfind(fld,' '))='_';
        info=setfield(info,fld,val);
    elseif(~isempty(strfind(line,'[')))
        %header skip
    end
    
end
fclose(fid);

if(isfield(info,'Wavelengths'))
    
    % fix a few strings that should be numeric
    info.Wavelengths=str2num(info.Wavelengths);
    info.Mod_Amp=str2num(info.Mod_Amp);
    info.Threshold=str2num(info.Threshold);
    info.ChanDis=str2num(info.ChanDis);
    
    % Fix the SD key
    if(strcmp(info.S_D_Key(end),',')); info.S_D_Key(end)=[]; end;
    keys=strsplit(info.S_D_Key,',');
    for idx=1:length(keys)
        a=strsplit(keys{idx},{'-',':'});
        info.SDkey(idx,1)=str2num(a{3});  % channel index
        info.SDkey(idx,2)=str2num(a{1});  % source
        info.SDkey(idx,3)=str2num(a{2});  % detector
    end
end

try
    % added to fix bug issue #57.  11/1/2018
    if isfield(info, 'ShortDetIndex') && isfield(info,'ShortBundles') && info.ShortBundles > 0
        info.ShortDetectors = length(str2num(info.ShortDetIndex));
    end
end

end