function snirf=NirsCoreData2SNIRF(data)

warning('off','MATLAB:fopen:FopenAllStagedRemoval');
snirf=SnirfClass;

for i=1:height(data.probe.link)
    ML(i)=MeasListClass();
    ML(i).detectorIndex=data.probe.link.detector(i);
    ML(i).sourceIndex=data.probe.link.source(i);
    ML(i).dataTypeLabel='raw';
    ML(i).dataType=1;
    ML(i).dataUnit='nm';
    ML(i).wavelengthActual=data.probe.link.type(i);
    ML(i).wavelengthIndex=find(ismember(data.probe.types,data.probe.link.type(i)));
end
  
DC = DataClass();
DC.time=data.time;
DC.dataTimeSeries=data.data;
DC.measurementList=ML;
snirf.data=DC;

probe=ProbeClass();
probe.wavelengths=data.probe.types;
probe.sourceLabels=data.probe.optodes(ismember(data.probe.optodes.Type,'Source'),:).Name;
probe.detectorLabels=data.probe.optodes(ismember(data.probe.optodes.Type,'Detector'),:).Name;
probe.landmarkLabels=data.probe.optodes(ismember(data.probe.optodes.Type,{'FID-anchor','FID-attractor'}),:).Name;

if(isa(data.probe,'nirs.core.Probe1020'))
    probe.sourcePos3D=[data.probe.optodes_registered(ismember(data.probe.optodes_registered.Type,'Source'),:).X...
        data.probe.optodes_registered(ismember(data.probe.optodes_registered.Type,'Source'),:).Y...
        data.probe.optodes_registered(ismember(data.probe.optodes_registered.Type,'Source'),:).Z];
    probe.detectorPos3D=[data.probe.optodes(ismember(data.probe.optodes_registered.Type,'Detector'),:).X...
        data.probe.optodes_registered(ismember(data.probe.optodes_registered.Type,'Detector'),:).Y...
        data.probe.optodes_registered(ismember(data.probe.optodes_registered.Type,'Detector'),:).Z];
    probe.landmarkPos3D=[data.probe.optodes(ismember(data.probe.optodes_registered.Type,{'FID-anchor','FID-attractor'}),:).X...
        data.probe.optodes_registered(ismember(data.probe.optodes_registered.Type,{'FID-anchor','FID-attractor'}),:).Y...
        data.probe.optodes_registered(ismember(data.probe.optodes_registered.Type,{'FID-anchor','FID-attractor'}),:).Z];

end

probe.sourcePos2D=[data.probe.optodes(ismember(data.probe.optodes.Type,'Source'),:).X...
    data.probe.optodes(ismember(data.probe.optodes.Type,'Source'),:).Y...
    data.probe.optodes(ismember(data.probe.optodes.Type,'Source'),:).Z];
probe.detectorPos2D=[data.probe.optodes(ismember(data.probe.optodes.Type,'Detector'),:).X...
    data.probe.optodes(ismember(data.probe.optodes.Type,'Detector'),:).Y...
    data.probe.optodes(ismember(data.probe.optodes.Type,'Detector'),:).Z];
probe.landmarkPos2D=[data.probe.optodes(ismember(data.probe.optodes.Type,{'FID-anchor','FID-attractor'}),:).X...
    data.probe.optodes(ismember(data.probe.optodes.Type,{'FID-anchor','FID-attractor'}),:).Y...
    data.probe.optodes(ismember(data.probe.optodes.Type,{'FID-anchor','FID-attractor'}),:).Z];

snirf.probe=probe;

for i=1:data.stimulus.count;
    stim=data.stimulus(data.stimulus.keys{i});
    snirf.stim(i)=StimClass();
    snirf.stim(i).name=data.stimulus.keys{i};
    snirf.stim(i).data=[stim.onset stim.dur stim.amp];
    snirf.stim(i).dataLabels={'Onset','Duration','Amplitude'};
    snirf.stim(i).states=[stim.onset stim.amp];
end

cnt=1;
for i=1:data.auxillary.count
    aux=data.auxillary(data.auxillary.keys{i});
    if(isa(aux,'nirs.core.GenericData'))
        snirf.aux(cnt)=AuxClass();
        snirf.aux(cnt).name=data.auxillary.keys{i};
        snirf.aux(cnt).dataTimeSeries=aux.data;
        snirf.aux(cnt).time=aux.time;
        cnt=cnt+1;
    end
end

if(~data.demographics.isempty)
    snirf.metaDataTags=MetaDataTagsClass();
    snirf.metaDataTags.tags=data.demographics.toStruct;
end