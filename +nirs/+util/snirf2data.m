function data=snirf2data(snirf)

dataTypeDefs{1,1}=1;
dataTypeDefs{1,2}='raw';
dataTypeDefs{2,1}=100;
dataTypeDefs{2,2}='raw-DC';
dataTypeDefs{3,1}=101;
dataTypeDefs{3,2}='raw-AC';
dataTypeDefs{4,1}=102;
dataTypeDefs{4,2}='raw-Phase';

data=nirs.core.Data;
data(:)=[];
cnt=1;
for i=1:length(snirf.nirs)
    tmpdata=nirs.core.Data;
    tmpdata(:)=[];
    for ii=1:length(snirf.nirs(i).data)
        tmpdata(ii).data=snirf.nirs(i).data(ii).dataTimeSeries;
        tmpdata(ii).time=snirf.nirs(i).data(ii).time;
        
        if(size(tmpdata(ii).time,2)>1)
            tmpdata(ii).time=tmpdata(ii).time';
        end
        if(size(tmpdata(ii).data,1)~=length(tmpdata(ii).time))
            tmpdata(ii).data=tmpdata(ii).data';
        end
        
        
        if(isfield(snirf.nirs(i).probe,'sourcePos3D'))
            tmpdata(ii).probe=nirs.core.Probe1020;
        else
            tmpdata(ii).probe=nirs.core.Probe;
        end
        
        if(~isfield(snirf.nirs(i).data(ii).measurementList(1),'dataTypeLabel'))
            for j=1:length(snirf.nirs(ii).data.measurementList)
                idx=find(ismember(vertcat(dataTypeDefs{:,1}),...
                    double(snirf.nirs(ii).data.measurementList(j).dataType)));
                if(isempty(idx))
                    warning(['unknown datatyoe: ' num2str(snirf.nirs(ii).data.measurementList(j).dataType)]);
                    idx=1;  % default to "raw"
                end
                snirf.nirs(i).data(ii).measurementList(j).dataTypeLabel=cellstr(dataTypeDefs{idx,2});
                
            end
        end
        
        if(~isfield(snirf.nirs(i).probe,'sourcePos') & isfield(snirf.nirs(i).probe,'sourcePos2D'))
            snirf.nirs(i).probe.sourcePos=snirf.nirs(i).probe.sourcePos2D;
        end
        if(~isfield(snirf.nirs(i).probe,'detectorPos') & isfield(snirf.nirs(i).probe,'detectorPos2D'))
            snirf.nirs(i).probe.detectorPos=snirf.nirs(i).probe.detectorPos2D;
        end
        
        if(~isfield(snirf.nirs(i).probe,'landmarkPos') & isfield(snirf.nirs(i).probe,'landmarkPos2D'))
            snirf.nirs(i).probe.landmarkPos=snirf.nirs(i).probe.landmarkPos2D;
        end
        
        if(size(snirf.nirs(i).probe.sourcePos,1)==2)
            snirf.nirs(i).probe.sourcePos=snirf.nirs(i).probe.sourcePos';
        end
        if(size(snirf.nirs(i).probe.detectorPos,1)==2)
            snirf.nirs(i).probe.detectorPos=snirf.nirs(i).probe.detectorPos';
            
        end
        
        if(~isfield(snirf.nirs(i).probe,'sourceLabels'))
            for j=1:size(snirf.nirs(i).probe.sourcePos,1)
                snirf.nirs(i).probe.sourceLabels{j,1}=['Source-' num2str(j)];
            end
        end
        
        if(~isfield(snirf.nirs(i).probe,'detectorLabels'))
            for j=1:size(snirf.nirs(i).probe.detectorPos,1)
                snirf.nirs(i).probe.detectorLabels{j,1}=['Detector-' num2str(j)];
            end
        end
        
        if(ischar(snirf.nirs(i).probe.sourceLabels))
            snirf.nirs(i).probe.sourceLabels={snirf.nirs(i).probe.sourceLabels};
        end
        if(ischar(snirf.nirs(i).probe.detectorLabels))
            snirf.nirs(i).probe.detectorLabels={snirf.nirs(i).probe.detectorLabels};
        end
        
        if(isfield(snirf.nirs(i).probe,'sourcePos3D') && size(snirf.nirs(i).probe.sourcePos3D,1)==3)
            snirf.nirs(i).probe.sourcePos3D=snirf.nirs(i).probe.sourcePos3D';
        end
        if(isfield(snirf.nirs(i).probe,'detectorPos3D') && size(snirf.nirs(i).probe.detectorPos3D,1)==3)
            snirf.nirs(i).probe.detectorPos3D=snirf.nirs(i).probe.detectorPos3D';
        end
        if(isfield(snirf.nirs(i).probe,'landmarkPos') && size(snirf.nirs(i).probe.landmarkPos,1)==2)
            snirf.nirs(i).probe.landmarkPos=snirf.nirs(i).probe.landmarkPos';
        end
        if(isfield(snirf.nirs(i).probe,'landmarkPos3D') && size(snirf.nirs(i).probe.landmarkPos3D,1)==3)
            snirf.nirs(i).probe.landmarkPos3D=snirf.nirs(i).probe.landmarkPos3D';
        end
        
        
        %make the optodes
        Type={}; X=[]; Y=[]; Z=[]; Units={}; Name={};
        for j=1:length(snirf.nirs(i).probe.sourceLabels)
            Type{end+1,1}='Source';
            Name{end+1,1}=snirf.nirs(i).probe.sourceLabels{j};
            Units{end+1,1}=snirf.nirs(i).metaDataTags.LengthUnit;
            X(end+1,1)=snirf.nirs(i).probe.sourcePos(j,1);
            Y(end+1,1)=snirf.nirs(i).probe.sourcePos(j,2);
            Z(end+1,1)=0;
        end
        for j=1:length(snirf.nirs(i).probe.detectorLabels)
            Type{end+1,1}='Detector';
            Name{end+1,1}=snirf.nirs(i).probe.detectorLabels{j};
            Units{end+1,1}=snirf.nirs(i).metaDataTags.LengthUnit;
            X(end+1,1)=snirf.nirs(i).probe.detectorPos(j,1);
            Y(end+1,1)=snirf.nirs(i).probe.detectorPos(j,2);
            Z(end+1,1)=0;
        end
        % add landmarks
        if(isfield(snirf.nirs(i).probe,'landmarkLabels') & isfield(snirf.nirs(i).probe,'landmarkPos'))
            for j=1:length(snirf.nirs(i).probe.landmarkLabels)
                Type{end+1,1}=snirf.nirs(i).probe.landmarkLabels{j}(strfind(snirf.nirs(i).probe.landmarkLabels{j},'FID'):end);
                Name{end+1,1}=snirf.nirs(i).probe.landmarkLabels{j}(1:strfind(snirf.nirs(i).probe.landmarkLabels{j},'FID')-1);
                Units{end+1,1}=snirf.nirs(i).metaDataTags.LengthUnit;
                X(end+1,1)=snirf.nirs(i).probe.landmarkPos(j,1);
                Y(end+1,1)=snirf.nirs(i).probe.landmarkPos(j,2);
                Z(end+1,1)=0;
            end
        end
        
        tmpdata(ii).probe.optodes=table(Name,X,Y,Z,Type,Units);
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
                if(size(snirf.nirs(i).probe.landmarkPos3D,2)>size(snirf.nirs(i).probe.landmarkPos3D,1));
                    snirf.nirs(i).probe.landmarkPos3D=snirf.nirs(i).probe.landmarkPos3D';
                end
                for j=1:length(snirf.nirs(i).probe.landmarkLabels)
                    Type{end+1,1}=snirf.nirs(i).probe.landmarkLabels{j}(strfind(snirf.nirs(i).probe.landmarkLabels{j},'FID'):end);
                    Name{end+1,1}=snirf.nirs(i).probe.landmarkLabels{j}(1:strfind(snirf.nirs(i).probe.landmarkLabels{j},'FID')-1);
                    Units{end+1,1}=snirf.nirs(i).metaDataTags.LengthUnit;
                    X(end+1,1)=snirf.nirs(i).probe.landmarkPos3D(j,1);
                    Y(end+1,1)=snirf.nirs(i).probe.landmarkPos3D(j,2);
                    Z(end+1,1)=snirf.nirs(i).probe.landmarkPos3D(j,3);
                end
            end
            tmpdata(ii).probe.optodes_registered=table(Name,X,Y,Z,Type,Units);
        end
        
        fds=fields(snirf.nirs(i).metaDataTags);
        for f=1:length(fds)
            a=snirf.nirs(i).metaDataTags.(fds{f});
            if(~isempty(a))
                tmpdata(ii).demographics(fds{f})=a;
            end
        end
        
        
        % now the link table
        source=[]; detector=[];
        if(isfield(snirf.nirs(i).probe,'wavelengths'))
            type=[];
        else
            type={};
        end
%<<<<<<< snirf-probe-distances
    end
    tmpdata(ii).probe.link=table(source,detector,type);
	% add registered optodes distances to probe fixeddistances - added by Peggy Skelly
	if(isfield(snirf.nirs(i).probe,'sourcePos3D'))
		tmpdata(ii).probe.fixeddistances=tmpdata(ii).probe.swap_reg.distances;
	end
	
    if(isfield(snirf.nirs(i),'stim'))
        for j=1:length(snirf.nirs(i).stim)
%=======
%        for j=1:length(snirf.nirs(i).data(ii).measurementList)
%            source(j,1)=snirf.nirs(i).data(ii).measurementList(j).sourceIndex;
%            detector(j,1)=snirf.nirs(i).data(ii).measurementList(j).detectorIndex;
%>>>>>>> master
            
            if(isfield(snirf.nirs(i).probe,'wavelengths') && ...
                    isfield(snirf.nirs(i).data(ii).measurementList(j),'wavelengthIndex'))
                type(j,1)=snirf.nirs(i).probe.wavelengths(snirf.nirs(i).data(ii).measurementList(j).wavelengthIndex);
            else
                type{j,1}=snirf.nirs(i).data(ii).measurementList(j).dataTypeLabel;
            end
        end
        tmpdata(ii).probe.link=table(source,detector,type);
        
        if(isfield(snirf.nirs(i),'stim'))
            for j=1:length(snirf.nirs(i).stim)
                
                if(size(snirf.nirs.stim(j).data,1)==length(snirf.nirs(i).stim(j).dataLabels))
                    snirf.nirs.stim(j).data=snirf.nirs.stim(j).data';
                end
                
                tbl=struct;
                for id=1:length(snirf.nirs(i).stim(j).dataLabels)
                    tbl=setfield(tbl,genvarname(snirf.nirs(i).stim(j).dataLabels{id}),snirf.nirs.stim(j).data(:,id));
                end
                
                st=nirs.design.StimulusEvents;
                st.name=snirf.nirs.stim(j).name;
                st.onset=snirf.nirs.stim(j).data(:,1);
                st.dur=snirf.nirs.stim(j).data(:,2);
                st.amp=snirf.nirs.stim(j).data(:,3);
                st.metadata=struct2table(tbl);
                tmpdata(ii).stimulus(st.name)=st;
            end
        end
        
        if(isfield(snirf.nirs(i),'aux'))
            for j=1:length(snirf.nirs(i).aux)
                st=nirs.core.GenericData;
                st.data=snirf.nirs(i).aux(j).dataTimeSeries;
                st.time=snirf.nirs(i).aux(j).time;
                
                tmpdata(ii).auxillary(snirf.nirs(i).aux(j).name)=st;
            end
        end
        if(isfield(snirf.nirs(i).probe,'frequencies'))
            tmpdata(ii).Fm=snirf.nirs(i).probe.frequencies(end);
        end
        
    end
    if(tmpdata(ii).Fm~=0)
        % deal with FD data
        data(i)=combinedata(tmpdata);
    else
        data=tmpdata(1);
    end
end



return


function data = combinedata(tmpdata)
data=tmpdata(1);




% phase=angle(data(i).data(:,lst));
%         data(i).data=abs(data(i).data);
%         data(i).data=[data(i).data phase];

return

