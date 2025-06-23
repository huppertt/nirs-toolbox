function data=readKernelMomentsSNIRF(filename)
% This function converts data read in using the Homer-3 SnirfClass to the
% nirs-toolbox nirs.core.Data structure



warning('off','MATLAB:fopen:FopenAllStagedRemoval')
snirfdata=SnirfLoad(filename);

data = nirs.core.Data;
data.description=filename;
data.demographics=Dictionary(snirfdata.metaDataTags.tags);
data.demographics('SNIRF Format')=snirfdata.formatVersion;


for stimIdx=1:length(snirfdata.stim)
    tbl=array2table(snirfdata.stim(stimIdx).data,'VariableNames',snirfdata.stim(stimIdx).dataLabels);

    stim=nirs.design.StimulusEvents;
    stim.name=snirfdata.stim(stimIdx).name;
    stim.metadata=tbl;
    stim.onset=tbl.Timestamp;
    stim.dur=tbl.Duration;
    stim.amp=tbl.Value;
    data.stimulus(stim.name)=stim;
end

for auxIdx=1:length(snirfdata.aux)
    aux = nirs.core.GenericData;
    aux.description=snirfdata.aux(auxIdx).name;
    aux.data=snirfdata.aux(auxIdx).dataTimeSeries;
    aux.time=snirfdata.aux(auxIdx).time;
    if(~isempty(snirfdata.aux(auxIdx).timeOffset))
        aux.time=aux.time+snirfdata.aux(auxIdx).timeOffset;
    end
    data.auxillary(snirfdata.aux(auxIdx).name)=aux;
end

data.auxillary('eeg') = getEEGfromAux(data);


data.time=snirfdata.data.time;
data.data=snirfdata.data.dataTimeSeries;

link=getML_link(snirfdata);

probe1020=nirs.core.Probe1020;
probe1020.link=link;

optodes.Name ={};
optodes.Type ={};
optodes.Units={};
for i=1:length(snirfdata.probe.sourceLabels)
    lab=['000' num2str(i)];
    lab=lab(end-3:end);
    optodes.Name{end+1,1} =['Source-' lab];
    optodes.Type{end+1,1}='Source';
    optodes.Units{end+1,1}='mm';
end
for i=1:length(snirfdata.probe.detectorLabels)
    lab=['000' num2str(i)];
    lab=lab(end-3:end);
    optodes.Name{end+1,1} =['Detector-' lab];
    optodes.Type{end+1,1}='Detector';
    optodes.Units{end+1,1}='mm';
end
optodes3D.Name=optodes.Name;
optodes3D.Type=optodes.Type;
optodes3D.Units=optodes.Units;
for i=1:length(snirfdata.probe.landmarkLabels)
    optodes3D.Name{end+1,1}=snirfdata.probe.landmarkLabels{i};
    optodes3D.Type{end+1,1}='FID-anchor';
    optodes3D.Units{end+1,1}='mm';
end



optodes.X = [   snirfdata.probe.sourcePos2D(:,1);...
    snirfdata.probe.detectorPos2D(:,1)];
optodes.Y = [   snirfdata.probe.sourcePos2D(:,2);...
    snirfdata.probe.detectorPos2D(:,2)];
optodes.Z = [   snirfdata.probe.sourcePos2D(:,3);...
    snirfdata.probe.detectorPos2D(:,3)];

optodes3D.X = [   snirfdata.probe.sourcePos3D(:,1);...
    snirfdata.probe.detectorPos3D(:,1);...
    snirfdata.probe.landmarkPos3D(:,1)];
optodes3D.Y = [   snirfdata.probe.sourcePos3D(:,2);...
    snirfdata.probe.detectorPos3D(:,2);...
    snirfdata.probe.landmarkPos3D(:,2)];
optodes3D.Z = [   snirfdata.probe.sourcePos3D(:,3);...
    snirfdata.probe.detectorPos3D(:,3);...
    snirfdata.probe.landmarkPos3D(:,3)];



for id=1:length(optodes.Type)
    idx=str2num(optodes.Name{id}(strfind(optodes.Name{id},'-')+1:end));
    if(strcmp(optodes.Type{id},'Source'))
        optodes.NameKernel{id,1}=lower(snirfdata.probe.sourceLabels{idx});
    elseif((strcmp(optodes.Type{id},'Detector')))
        optodes.NameKernel{id,1}=lower(snirfdata.probe.detectorLabels{idx});
    else
        optodes.NameKernel{id,1}=optodes.Name{id};
    end
end
for id=1:length(optodes3D.Type)
    idx=str2num(optodes3D.Name{id}(strfind(optodes3D.Name{id},'-')+1:end));
    if(strcmp(optodes3D.Type{id},'Source'))
        optodes3D.NameKernel{id,1}=lower(snirfdata.probe.sourceLabels{idx});
    elseif((strcmp(optodes3D.Type{id},'Detector')))
        optodes3D.NameKernel{id,1}=lower(snirfdata.probe.detectorLabels{idx});
    else
        optodes3D.NameKernel{id,1}=optodes3D.Name{id};
    end
end

probe1020.optodes_registered=struct2table(optodes3D);
probe1020.optodes=struct2table(optodes);

data.probe=probe1020;

for idx=1:height(data.probe.link); 
    NN{idx,1}=strcat(snirfdata.probe.sourceLabels{data.probe.link.source(idx)},...
        '_',snirfdata.probe.detectorLabels{data.probe.link.detector(idx)}); 
end;
data.probe.link.NameKernel=lower(NN);

function link=getML_link(snirfdata)

warning('off','MATLAB:structOnObject')
s=struct(snirfdata.data.measurementList(1));
for i=1:length(snirfdata.data.measurementList);
    s(i)=struct(snirfdata.data.measurementList(i));
end;
link=struct2table(s);

% remove empty
flds=link.Properties.VariableNames;
for i=1:length(flds)
    try
        if(isempty(vertcat(link.(flds{i}){:})))
            link.(flds{i})=[];
        end
    end
    try
        if(isempty(vertcat(link.(flds{i})(:))))
            link.(flds{i})=[];
        end
    end
end
link.supportedFomats=[];
link.err=[];

type=[];
for i=1:height(link);
    if(link.dataType==301);
        dataType=301;
        type(i,1)=snirfdata.probe.wavelengths(link.wavelengthIndex(i));
        moment(i,1)=link.dataTypeIndex(i);
    elseif(link.dataType==99999);
        dataType=99999;
        type{i,1}=link.dataTypeLabel{i};
    end

end

link2.source=link.sourceIndex;
link2.detector=link.detectorIndex;
link2.type=type;
if(dataType==301)
    link2.moment=moment;
end
link=struct2table(link2);




function eegdata = getEEGfromAux(data)

keys=data.auxillary.keys';
eegnames={};
for i=1:length(keys)
    if(contains(keys{i},'eeg') & contains(keys{i},'microV'))
        eegnames{end+1,1}=keys{i};
    end
end
eegpos=strrep(eegnames,'eeg_','');
eegpos=strrep(eegpos,'_microV','');

for i=1:length(eegnames)
    aux=data.auxillary(eegnames{i});
    d(:,i)=aux.data;
end

eegdata=eeg.core.Data;
eegdata.time=aux.time;
eegdata.data=d;
eegdata.description=data.description;
eegdata.stimulus=data.stimulus;
eegdata.demographics=data.demographics;

eegdata.probe=eeg.core.Probe(eegpos);

return

