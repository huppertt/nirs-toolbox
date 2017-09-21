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
if(isempty(file)); raw=[]; return; end;
info = parsehdr(fullfile(folder,file(1).name));

%now read the probeInfo file and convert to a probe class
file = dir(fullfile(folder,'*_probeInfo.mat'));
if(isempty(file)); raw=[]; return; end;

load(fullfile(folder,file(1).name)); % --> probeInfo

SrcPos = probeInfo.probes.coords_s2;
SrcPos(:,3)=0;

DetPos = probeInfo.probes.coords_d2;
DetPos(:,3)=0;

[s,d]=find(info.S_D_Mask);

%% Fix added based on bug from Guilherme A. Zimeo Morais created an issue 2016-08-24
[s,a]=sort(s); d=d(a);

link=table(repmat(s,length(info.Wavelengths),1),...
    repmat(d,length(info.Wavelengths),1),...
    reshape(repmat(info.Wavelengths(:)',length(s),1),[],1),...
    'VariableNames',{'source','detector','type'});

probe = nirs.core.Probe(SrcPos,DetPos,link);

% Not sure why the units on the 2D probe in the NIRx file are so off
scale=mean(info.ChanDis(:)./probe.distances(1:length(info.ChanDis(:))));
probe.optodes.X=scale*probe.optodes.X;
probe.optodes.Y=scale*probe.optodes.Y;
probe.optodes.Z=scale*probe.optodes.Z;

%% Now, let's get the data
raw = nirs.core.Data();





lst=find(ismember(info.SDkey(:,2:3),[probe.link.source probe.link.detector],'rows'));
for idx=1:length(info.Wavelengths)
    file = dir(fullfile(folder,['*.wl' num2str(idx)]));
    
    if(file(1).bytes==0)
        error('zero byte file');
    end
    
    d = dlmread(fullfile(folder,file(1).name));
    raw.data=[raw.data d(:,lst)];
end

raw.time=[0:size(raw.data,1)-1]/info.SamplingRate;

raw.description=info.FileName;

% Add the demographics info
file = dir(fullfile(folder,'*.inf'));
demoinfo = parsehdr(fullfile(folder,file(1).name));
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
if(length(info.ChanDis)==length(s)*2)
    raw.demographics('hyperscan')=info.FileName;
end


% Now add stimulus info
if(~isempty(info.Events))
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
fidS=table(fid_1020.Name(probeInfo.probes.index_s(:,2)),...
    probe.srcPos(probeInfo.probes.index_s(:,1),1),...
    probe.srcPos(probeInfo.probes.index_s(:,1),2),...
    probe.srcPos(probeInfo.probes.index_s(:,1),3),...
    repmat({'FID-anchor'},length(probeInfo.probes.index_s),1),repmat({'mm'},length(probeInfo.probes.index_s),1),...
    'VariableNames',{'Name','X','Y','Z','Type','Units'});
fidD=table(fid_1020.Name(probeInfo.probes.index_d(:,1)),...
    probe.detPos(probeInfo.probes.index_d(:,1),1),...
    probe.detPos(probeInfo.probes.index_d(:,1),2),...
    probe.detPos(probeInfo.probes.index_d(:,1),3),...
    repmat({'FID-anchor'},length(probeInfo.probes.index_d),1),repmat({'mm'},length(probeInfo.probes.index_d),1),...
    'VariableNames',{'Name','X','Y','Z','Type','Units'});

% and concatinate it to the probe
probe.optodes=[probe.optodes; fidS; fidD];

if(registerprobe)
    probe1020=nirs.util.registerprobe1020(probe);
    
    
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
    
    % This will allow NIRFAST to directly use the info for the BEM model
    lambda=unique(probe.link.type);
    prop{1} = nirs.media.tissues.skin(lambda);
    prop{2} = nirs.media.tissues.bone(lambda);
    prop{3} = nirs.media.tissues.water(lambda);
    prop{4} = nirs.media.tissues.brain(0.7, 50,lambda);
    
    fwdBEM=nirs.forward.NirfastBEM;
    fwdBEM.mesh=BEM;
    fwdBEM.prop  = prop;
    
    probe1020=probe1020.regsister_mesh2probe(fwdBEM.mesh);
    probe1020.opticalproperties=prop;
    
    
    SrcPos3D = probeInfo.probes.coords_s3;
    DetPos3D = probeInfo.probes.coords_d3;
    FID3D = [probeInfo.geom.NIRxHead.ext1020sys.coords3d(probeInfo.probes.index_s(:,2),:);...
        probeInfo.geom.NIRxHead.ext1020sys.coords3d(probeInfo.probes.index_d(:,2),:)];
    
    XYZ=[SrcPos3D; DetPos3D; FID3D];
    
    fidPts=probe1020.optodes_registered(ismember(probe1020.optodes_registered.Type,'FID-anchor'),:);
    XYZ_reg=[fidPts.X fidPts.Y fidPts.Z];
    XYZ_reg(:,4)=1;
    FID3D(:,4)=1;
    XYZ(:,4)=1;
    R=FID3D\XYZ_reg;
    
    XYZ=XYZ*R;
    
    probe1020.optodes_registered=probe1020.optodes;
    probe1020.optodes_registered.X=XYZ(:,1);
    probe1020.optodes_registered.Y=XYZ(:,2);
    probe1020.optodes_registered.Z=XYZ(:,3);
    
    raw.probe=probe1020;
else
    raw.probe=probe;
end

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
        if(isempty(strfind(line,'#')))
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
                if(isempty(strfind(line,'#')))
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

end


