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

disp(['Loading : ' folder]);
probe=nirs.io.loadNIRxProbe(folder,registerprobe);


% this is a hack that allowed me to fix a study where the NIRx aquistion had
% been setup wrong and was missing some channels.  NIRx saves all data in
% the wl files (even those not part of the probe).  Creating a SDmask.mat
% file overwrites the info in the NIRx file
if(~isempty(dir(fullfile(folder,'*SDmask.mat'))))
    fi=dir(fullfile(folder,'*SDmask.mat'));
    load(fullfile(folder,fi(1).name));
    disp('loading mask from file');
    info.S_D_Mask=SD_mask;
    
    
    [s,d]=find(info.S_D_Mask);
    
    
    %% Fix added based on bug from Guilherme A. Zimeo Morais created an issue 2016-08-24
    [s,a]=sort(s); d=d(a);
    
    link=table(s,d,...
        'VariableNames',{'source','detector'});
    probe.link=link;
    
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


if(isfield(info,'ShortDetectors') & info.Detectors>size(probe.detPos,1))
    % this is an intermediate version of the software when short
    % seperations were first introduced.
      if (info.ShortDetectors >= 1)
        lst=find(ismember(probe.optodes.Type,{'Source'}));
        lst2=find(ismember(probe.optodes.Type,{'Detector'}));
        lst3=find(~ismember(probe.optodes.Type,{'Source','Detector'}));
    
        DetShort =probe.optodes(lst,:);
        DetShort3D =probe.optodes_registered(lst,:);
        for i=1:length(lst)
            str=['000' num2str(i+info.Detectors-info.ShortDetectors)];
            DetShort.Name{i}=['Detector-' str(end-3:end)];
            DetShort.Type{i}='Detector';
            DetShort.X(i)=DetShort.X(i)+eps(1);
            DetShort3D.Name{i}=['Detector-' str(end-3:end)];
            DetShort3D.Type{i}='Detector';
            DetShort3D.X(i)=DetShort3D.X(i)+1;
            
        end
        probe.optodes=[probe.optodes(lst2(1:info.Detectors-info.ShortDetectors),:); ...
            DetShort; probe.optodes(lst,:); probe.optodes(lst3,:)];
        probe.optodes_registered=[probe.optodes_registered(lst2(1:info.Detectors-info.ShortDetectors),:); ...
            DetShort3D; probe.optodes_registered(lst,:); probe.optodes_registered(lst3,:)];
        
        info.S_D_Mask(:,info.Detectors-info.ShortDetectors+1:end)=eye(info.ShortDetectors);
      end
        
end



%% Now, let's get the data
file = rdir(fullfile(folder,'*.nirs' ));
if(~isempty(file))
    raw=nirs.io.loadDotNirs(file(1).name,true);
    probe.link=raw.probe.link;
    raw.probe=probe;
    
else
    [s,d]=find(info.S_D_Mask);
    [s,a]=sort(s); d=d(a);
    
    link=table(s,d,...
        'VariableNames',{'source','detector'});
    if(height(probe.link)<height(link))
        probe.link=link;
    end
    
    if(info.Wavelengths==2)
        % not sure why
        info.Wavelengths=[780 850];
    end
    
    probe.link=[[probe.link; probe.link] ...
        table(reshape(repmat(info.Wavelengths(:)',...
        height(probe.link),1),[],1),'VariableNames',{'type'})];
    
    if(~isfield(info,'SDkey') & isfield(info,'S_D_Key'))
        info.SDkey=info.S_D_Key;
    end
    

    
    
    % read the standard wl files
    raw = nirs.core.Data();
    raw.probe=probe;
    
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

raw.probe.fixeddistances=raw.probe.swap_reg.distances;


if(isfield(info,'ShortDetectors') && ~isempty(info.ShortDetectors) && ...
        info.ShortDetectors) 
    
   j=nirs.modules.LabelShortSeperation;
   j.max_distance=10;
   raw=j.run(raw);
   % raw.probe.link(ismember(raw.probe.link.detector,info.ShortDetIndex),:)

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
if(exist('demoinfo'))
flds=fields(demoinfo);
for idx=1:length(flds)
    demo(flds{idx})=demoinfo.(flds{idx});
end
end
raw.demographics=demo;



% There is a slight error in the NIRx files for hyperscanning that I am
% going to exploit to add this info.  For hyperscanning files, the
% info.S_D_Mask covers only 1 subject but the info.SD_Key field is the full
% (both subjects) length.
if(isfield(info,'ChanDis'))
    if(length(info.ChanDis)==height(probe.link))
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

if(isfield(info,'Wavelengths') & isfield(info,'ChanDis'))
    
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