function [snirf,snirfarray] = validateSNIRF(data)

snirf.formatVersion='1.0';

for i=1:length(data)
    
    snirf.nirs(i).metaDataTags.SubjectID=nirs.util.dictionaryGet(data(i).demographics,{'subjid','subject','name'},'annonymous');
    snirf.nirs(i).metaDataTags.MeasurementDate=nirs.util.dictionaryGet(data(i).demographics,{'MeasurementDate'},'?');
    snirf.nirs(i).metaDataTags.MeasurementTime=nirs.util.dictionaryGet(data(i).demographics,{'MeasurementTime'},'?');
    snirf.nirs(i).metaDataTags.LengthUnit=data(i).probe.optodes.Units{1};
    snirf.nirs(i).metaDataTags.TimeUnit='s';
    
    %     snirf.nirs(i).metaDataTags.AccessionNumber =
    %     snirf.nirs(i).metaDataTags.InstanceNumber =
    %     snirf.nirs(i).metaDataTags.UnixTime =
    
    for j=1:data(i).demographics.count
        snirf.nirs(i).metaDataTags=setfield(snirf.nirs(i).metaDataTags,data(i).demographics.keys{j},data(i).demographics(data(i).demographics.keys{j}));
    end
    
    snirf.nirs(i).metaDataTags.filedescription=data(i).description;
    snirf.nirs(i).metaDataTags.SNIRF_createDate=datestr(now,'yyyy-mm-dd');
    snirf.nirs(i).metaDataTags.SNIRF_createTime=[datestr(now,'hh:mm:ss') 'Z'];
    
    if(data(i).Fm==0)
        %CW-NIRS
        snirf.nirs(i).data(1).dataTimeSeries=data(i).data';
        snirf.nirs(i).data(1).time=data(i).time;
        
        wavelengths=unique(data(i).probe.link.type);
        
        for k=1:height(data(i).probe.link);
            snirf.nirs(i).data(1).measurementList(k).sourceIndex=data(i).probe.link.source(k);
            snirf.nirs(i).data(1).measurementList(k).detectorIndex=data(i).probe.link.detector(k);
            snirf.nirs(i).data(1).measurementList(k).wavelengthIndex=find(ismember(wavelengths,data(i).probe.link.type(k)));
            
            snirf.nirs(i).data(1).measurementList(k).dataType = 1;
            snirf.nirs(i).data(1).measurementList(k).dataTypeIndex = 1;
            if(isnumeric(data(i).probe.link.type(k)))
                snirf.nirs(i).data(1).measurementList(k).dataTypeLabel='raw';
            else
                snirf.nirs(i).data(1).measurementList(k).dataTypeLabel=data(i).probe.link.type(k);
                if(strcmp(lower(snirf.nirs(i).data(1).measurementList(k).dataTypeLabel),'hbo'))
                    snirf.nirs(i).data(1).measurementList(k).dataTypeLabel='HbO';
                end
                if(strcmp(lower(snirf.nirs(i).data(1).measurementList(k).dataTypeLabel),'hbr'))
                    snirf.nirs(i).data(1).measurementList(k).dataTypeLabel='HbR';
                end
            end
            
            %        snirf.nirs(i).data(1).measurementList(k).sourcePower
            %        snirf.nirs(i).data(1).measurementList(k).detectorGain
            %       snirf.nirs(i).data(1).measurementList(k).moduleIndex
        end
    else
        % FD-NIRS
        ModFreq=unique(data(i).probe.link.ModFreq);
        lst=find(data(i).probe.link.ModFreq~=0);
        phase=angle(data(i).data(:,lst));
        data(i).data=abs(data(i).data);
        data(i).data=[data(i).data phase];
        l=data(i).probe.link(lst,:);
        l.ModFreq=-l.ModFreq;
        data(i).probe.link=[data(i).probe.link; l];
        
        
        
        wavelengths=unique(data(i).probe.link.type);
        
        ii=1;
        snirf.nirs(i).data(ii).time=data(i).time;
        link=data(i).probe.link;
        snirf.nirs(i).data(ii).dataTimeSeries=data(i).data';
        
        
        for k=1:height(link);
            snirf.nirs(i).data(ii).measurementList(k).sourceIndex=link.source(k);
            snirf.nirs(i).data(ii).measurementList(k).detectorIndex=link.detector(k);
            snirf.nirs(i).data(ii).measurementList(k).wavelengthIndex=find(ismember(wavelengths,link.type(k)));
            
            if(data(i).probe.link.ModFreq(k)==0)
                snirf.nirs(i).data(ii).measurementList(k).dataType = 1;  % DC-CW
                snirf.nirs(i).data(ii).measurementList(k).dataTypeLabel='raw-DC';
            elseif(data(i).probe.link.ModFreq(k)>0)
                snirf.nirs(i).data(ii).measurementList(k).dataType = 101;  % FD-Amp
                snirf.nirs(i).data(ii).measurementList(k).dataTypeLabel='raw-AC';
            else
                snirf.nirs(i).data(ii).measurementList(k).dataType = 102;  % FD-Phase
                snirf.nirs(i).data(ii).measurementList(k).dataTypeLabel='raw-Phase';
            end
            snirf.nirs(i).data(ii).measurementList(k).dataTypeIndex = 1;
            
            
            %        snirf.nirs(i).data(1).measurementList(k).sourcePower
            %        snirf.nirs(i).data(1).measurementList(k).detectorGain
            snirf.nirs(i).data(ii).measurementList(k).moduleIndex=find(ismember(ModFreq,abs(link.ModFreq(k))));
        end
        snirf.nirs(i).probe.frequencies=ModFreq;
        
        
    end
    
    cnt=1;
    for k=1:data(i).stimulus.count
        st=data(i).stimulus(data(i).stimulus.keys{k});
        if(isa(st,'nirs.design.StimulusVector'))
            a=nirs.core.GenericData;
            a.data=st.vector;
            a.time=st.time;
            data(i).auxillary(data(i).stimulus.keys{k})=a;
        else
            snirf.nirs(i).stim(cnt).name=st.name;
            snirf.nirs(i).stim(cnt).data =[st.onset st.dur st.amp]';
            cnt=cnt+1;
        end
    end
    
    if(isnumeric(wavelengths))
        snirf.nirs(i).probe.wavelengths =wavelengths;
    end
    %snirf.nirs(i).probe.wavelengthsEmission
    
    
    snirf.nirs(i).probe.sourcePos2D = data(i).probe.srcPos(:,1:2)';
    snirf.nirs(i).probe.sourceLabels = {data(i).probe.optodes(ismember(data(i).probe.optodes.Type,'Source'),:).Name{:}}';
    snirf.nirs(i).probe.detectorPos2D = data(i).probe.detPos(:,1:2)';
    snirf.nirs(i).probe.detectorLabels = {data(i).probe.optodes(ismember(data(i).probe.optodes.Type,'Detector'),:).Name{:}}';
    
    if(isa( data(i).probe,'nirs.core.Probe1020'))
        p=data(i).probe.swap_reg;
        snirf.nirs(i).probe.sourcePos3D=p.srcPos';
        snirf.nirs(i).probe.detectorPos3D=p.detPos';
    end
    
    
    %    snirf.nirs(i).probe.timeDelays
    %    snirf.nirs(i).probe.timeDelayWidths
    %    snirf.nirs(i).probe.momentOrders
    %    snirf.nirs(i).probe.correlationTimeDelays
    %    snirf.nirs(i).probe.correlationTimeDelayWidths
    
    lst=find(ismember(data(i).probe.optodes.Type,{'FID-anchor','FID-attractor'}));
    if(~isempty(lst))
        snirf.nirs(i).probe.landmarkPos2D=[data(i).probe.optodes.X(lst) data(i).probe.optodes.Y(lst)]';
        if(isa( data(i).probe,'nirs.core.Probe1020'))
            snirf.nirs(i).probe.landmarkPos3D=[data(i).probe.optodes_registered.X(lst) ...
                data(i).probe.optodes_registered.Y(lst) data(i).probe.optodes_registered.Z(lst)]';
            
        end
        snirf.nirs(i).probe.landmarkLabels=strcat({data.probe.optodes.Name{lst}},...
            repmat({' '},1,length(lst)), {data.probe.optodes.Type{lst}})';
    end
    
    snirf.nirs(i).probe.useLocalIndex=1;
    
    
    if(data(i).auxillary.count>0)
        cnt=1;
        for j=1:data(i).auxillary.count
            st=data(i).auxillary(data(i).auxillary.keys{j});
            if(isa(st,'nirs.core.GenericData'))
                snirf.nirs(i).aux(cnt).name=data(i).auxillary.keys{j};
                snirf.nirs(i).aux(cnt).dataTimeSeries=st.data;
                snirf.nirs(i).aux(cnt).time=st.time;
                snirf.nirs(i).aux(cnt).timeOffset=0;
                cnt=cnt+1;
            end
        end
    end
    
    
    
    
end

if(nargout==2)
    array=makearray(snirf);
    flds=fields(array);
    snirfarray = cell(length(flds),2);
    for j=1:length(flds)
        str=strsplit(flds{j},'__');
        n='';
        for k=2:length(str)
            n=[n '/' str{k}];
        end
        snirfarray{j,1}=n;
        snirfarray{j,2}=array.(flds{j});
    end
    
end
end

function array = makearray(snirf,array,str)

if(isfield(snirf,'nirs'))
    
    for i=1:length(snirf.nirs)
        if(length(snirf.nirs(i).data)==1)
            snirf.nirs(i).data1=snirf.nirs(i).data;
            a(i)=rmfield(snirf.nirs(i),'data');
        else
            a=snirf.nirs;
        end
    end
    snirf.nirs=a;
    clear a;
end

if(isfield(snirf,'nirs'))
    for i=1:length(snirf.nirs)
        if(isfield(snirf.nirs(i),'aux'))
            if(length(snirf.nirs(i).aux)==1)
                snirf.nirs(i).aux1=snirf.nirs(i).aux;
                ab(i)=rmfield(snirf.nirs(i),'aux');
            else
                ab(i)=snirf.nirs(i);
            end
        else
            ab(i)=snirf.nirs(i);
        end
    end
    snirf.nirs=ab;
    clear ab;
end


if(isfield(snirf,'stim'))
    if(length(snirf.stim)==1)
        snirf.stim1=snirf.stim;
        snirf=rmfield(snirf,'stim');
    end
end

if(nargin<2)
    array=struct;
end
if(nargin<3)
    str='snirf';
end

for i=1:length(snirf)
    flds = fields(snirf(i));
    
    if(length(snirf)==1)
        s2=[];
    else
        s2=num2str(i);
    end
    
    for j=1:length(flds)
        if(isstruct(snirf(i).(flds{j})))
            a = makearray(snirf(i).(flds{j}),array,[str s2 '__' flds{j} ]);
            ff=fields(a);
            for k=1:length(ff)
                array = setfield(array,ff{k},a.(ff{k}));
            end
            
            
        elseif((iscell(snirf(i).(flds{j})) & length(snirf(i).(flds{j}))>1));
            if(~isempty(strfind(flds{j},'Labels')))
                 array = setfield(array,[str s2 '__' flds{j}],snirf(i).(flds{j}));
            else
            for k=1:length(snirf(i).(flds{j}))
                array = setfield(array,[str s2 '__' flds{j} num2str(k)],snirf(i).(flds{j}){k});
            end
            end
        else
            array = setfield(array,[str s2 '__' flds{j}],snirf(i).(flds{j}));
        end
    end
end
end



% /nirs(i)/metaDataTags/SubjectID - string
% /nirs(i)/metaDataTags/MeasurementDate (YYYY-MM-DD) . string
% /nirs(i)/metaDataTags/MeasurementTime (hh:mm:ss.sTZD) string
% /nirs(i)/metaDataTags/LengthUnit (e.g. mm,um,cm) . string
% /nirs(i)/metaDataTags/TimeUnit (e.g. s, us).  string

%"AccessionNumber" is similar to the DICOM tag "Accession Number"(0008,0050), as defined in the DICOM standard (ISO 12052)
%"InstanceNumber" is defined similarly to the DICOM tag "Instance Number" (0020,0013
%"UnixTime" defines the Unix Epoch Time
%
%
%
% Supported measurementList(k).dataType values in dataTimeSeries
% 001-100: Raw - Continuous Wave (CW)
%
% 001 - Amplitude
% 051 - Fluorescence Amplitude
% 101-200: Raw - Frequency Domain (FD)
%
% 101 - AC Amplitude
% 102 - Phase
% 151 - Fluorescence Amplitude
% 152 - Fluorescence Phase
% 201-300: Raw - Time Domain - Gated (TD Gated)
%
% 201 - Amplitude
% 251 - Fluorescence Amplitude
% 301-400: Raw - Time domain ? Moments (TD Moments)
%
% 301 - Amplitude
% 351 - Fluorescence Amplitude
% 401-500: Raw - Diffuse Correlation Spectroscopy (DCS):
%
% 401 - g2
% 410 - BFi
% 99999: Processed
%
% Supported measurementList(k).dataTypeLabel values in dataTimeSeries
% Tag Name	Meanings
% "dOD"	Change in optical density
% "mua"	Absorption coefficient
% "musp"	Scattering coefficient
% "HbO"	Oxygenated hemoglobin (oxyhemoglobin) concentration
% "HbR"	Deoxygenated hemoglobin (deoxyhemoglobin) concentration
% "HbT"	Total hemoglobin concentration
% "H2O"	Water content
% "Lipid"	Lipid concentration
% "BFi"	Blood flow index
% "HRF dOD"	Hemodynamic response function for change in optical density
% "HRF HbO"	Hemodynamic response function for oxyhemoglobin concentration
% "HRF HbR"	Hemodynamic response function for deoxyhemoglobin concentration
% "HRF HbT"	Hemodynamic response function for total hemoglobin concentration
% "HRF BFi"	Hemodynamic response function for blood flow index