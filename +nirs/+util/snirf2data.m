function data=snirf2data(snirf)

data=nirs.core.Data;
data(:)=[];
for i=1:length(snirf.nirs)
    data(i)=nirs.core.Data;
    data(i).data=snirf.nirs(i).data.dataTimeSeries;
    data(i).time=snirf.nirs(i).data.time;
    
    if(isfield(snirf.nirs(i).probe,'sourcePos3D'))
        data(i).probe=nirs.core.Probe1020;
    else
        data(i).probe=nirs.core.Probe;
    end
    
    if(ischar(snirf.nirs.probe.sourceLabels))
        snirf.nirs.probe.sourceLabels={snirf.nirs.probe.sourceLabels};
    end
    if(ischar(snirf.nirs.probe.detectorLabels))
        snirf.nirs.probe.detectorLabels={snirf.nirs.probe.detectorLabels};
    end
    
    
    %make the optodes
    Type={}; X=[]; Y=[]; Z=[]; Units={}; Name={};
    for j=1:length(snirf.nirs.probe.sourceLabels)
        Type{end+1,1}='Source';
        Name{end+1,1}=snirf.nirs(i).probe.sourceLabels{j};
        Units{end+1,1}=snirf.nirs(i).metaDataTags.LengthUnit;
        X(end+1,1)=snirf.nirs(i).probe.sourcePos(j,1);
        Y(end+1,1)=snirf.nirs(i).probe.sourcePos(j,2);
        Z(end+1,1)=snirf.nirs(i).probe.sourcePos(j,3);
    end
    for j=1:length(snirf.nirs(i).probe.detectorLabels)
        Type{end+1,1}='Detector';
        Name{end+1,1}=snirf.nirs(i).probe.detectorLabels{j};
        Units{end+1,1}=snirf.nirs(i).metaDataTags.LengthUnit;
        X(end+1,1)=snirf.nirs(i).probe.detectorPos(j,1);
        Y(end+1,1)=snirf.nirs(i).probe.detectorPos(j,2);
        Z(end+1,1)=snirf.nirs(i).probe.detectorPos(j,3);
    end
    % add landmarks
    if(isfield(snirf.nirs(i).probe,'landmarkLabels'));
        for j=1:length(snirf.nirs(i).probe.landmarkLabels)
            Type{end+1,1}=snirf.nirs(i).probe.landmarkLabels{j}(strfind(snirf.nirs(i).probe.landmarkLabels{j},'FID'):end);
            Name{end+1,1}=snirf.nirs(i).probe.landmarkLabels{j}(1:strfind(snirf.nirs(i).probe.landmarkLabels{j},'FID')-1);
            Units{end+1,1}=snirf.nirs(i).metaDataTags.LengthUnit;
            X(end+1,1)=snirf.nirs(i).probe.landmarkPos(j,1);
            Y(end+1,1)=snirf.nirs(i).probe.landmarkPos(j,2);
            Z(end+1,1)=snirf.nirs(i).probe.landmarkPos(j,3);
        end
    end
    
    data(i).probe.optodes=table(Name,X,Y,Z,Type,Units);
    if(isfield(snirf.nirs(i).probe,'sourcePos3D'))
        %make the optodes_registered
        Type={}; X=[]; Y=[]; Z=[]; Units={}; Name={};
        for j=1:length(snirf.nirs.probe.sourceLabels)
            Type{end+1,1}='Source';
            Name{end+1,1}=snirf.nirs(i).probe.sourceLabels{j};
            Units{end+1,1}=snirf.nirs(i).metaDataTags.LengthUnit;
            X(end+1,1)=snirf.nirs(i).probe.sourcePos3D(j,1);
            Y(end+1,1)=snirf.nirs(i).probe.sourcePos3D(j,2);
            Z(end+1,1)=snirf.nirs(i).probe.sourcePos3D(j,3);
        end
        for j=1:length(snirf.nirs(i).probe.detectorLabels)
            Type{end+1,1}='Detector';
            Name{end+1,1}=snirf.nirs(i).probe.detectorLabels{j};
            Units{end+1,1}=snirf.nirs(i).metaDataTags.LengthUnit;
            X(end+1,1)=snirf.nirs(i).probe.detectorPos3D(j,1);
            Y(end+1,1)=snirf.nirs(i).probe.detectorPos3D(j,2);
            Z(end+1,1)=snirf.nirs(i).probe.detectorPos3D(j,3);
        end
        if(isfield(snirf.nirs(i).probe,'landmarkLabels'));
            % add landmarks
            for j=1:length(snirf.nirs(i).probe.landmarkLabels)
                Type{end+1,1}=snirf.nirs(i).probe.landmarkLabels{j}(strfind(snirf.nirs(i).probe.landmarkLabels{j},'FID'):end);
                Name{end+1,1}=snirf.nirs(i).probe.landmarkLabels{j}(1:strfind(snirf.nirs(i).probe.landmarkLabels{j},'FID')-1);
                Units{end+1,1}=snirf.nirs(i).metaDataTags.LengthUnit;
                X(end+1,1)=snirf.nirs(i).probe.landmarkPos3D(j,1);
                Y(end+1,1)=snirf.nirs(i).probe.landmarkPos3D(j,2);
                Z(end+1,1)=snirf.nirs(i).probe.landmarkPos3D(j,3);
            end
        end
        data(i).probe.optodes_registered=table(Name,X,Y,Z,Type,Units);
    end
    
    fds=fields(snirf.nirs(i).metaDataTags);
    for f=1:length(fds)
        data(i).demographics(fds{f})=snirf.nirs(i).metaDataTags.(fds{f});
    end
    
    
    % now the link table
    source=[]; detector=[]; type=[];
    for j=1:length(snirf.nirs(i).data.measurementList)
        source(j,1)=snirf.nirs.data.measurementList(j).sourceIndex;
        detector(j,1)=snirf.nirs.data.measurementList(j).detectorIndex;
        type(j,1)=snirf.nirs.probe.wavelengths(snirf.nirs.data.measurementList(j).wavelengthIndex);
    end
    data(i).probe.link=table(source,detector,type);
    
    for j=1:length(snirf.nirs(i).stim)
        st=nirs.design.StimulusEvents;
        st.name=snirf.nirs.stim(j).name;
        st.onset=snirf.nirs.stim(j).data(:,1);
        st.dur=snirf.nirs.stim(j).data(:,2);
        st.amp=snirf.nirs.stim(j).data(:,3);
        data(i).stimulus(st.name)=st;
    end
    
    
    
end



return
