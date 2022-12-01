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

        %Used to make the optodes tables
        Type={}; X=[]; Y=[]; Z=[]; Units={}; Name={};

        % Get best matches for source, detector, landmark positions
        %    rotate if needed, detect ambiguous 2/3 probe detector
        %    configurations (if possible)
        posTypes={'source','detector','landmark'};

        for p = 1:length(posTypes)
            posType=posTypes{p};
            posLabelField=strcat(posType,'Labels');
            posOrigLabelField=strcat(posType,'OrigLabels');
            posField=strcat(posType,'Pos');
            pos2DField=strcat(posType,'Pos2D');
            pos3DField=strcat(posType,'Pos3D');
            nField=strcat(posType,'N');
            posCapital=strcat(upper(posType(1)),posType(2:end));
            n=nan;
            % Mark if fields need to be transposed
            rotateFields=nan;
            %Mark if position field was assigned
            posFieldPreexists=~isfield(snirf.nirs(i).probe,posField);

            % number of elements most clearly defined by number of labels
            if(isfield(snirf.nirs(i).probe,posLabelField))
                n=max(size(snirf.nirs(i).probe.(posLabelField)));

                 % If label exists but is not a cell string, convert it
                if(ischar(snirf.nirs(i).probe.(posLabelField)))
                    snirf.nirs(i).probe.(posLabelField)=cellstr(snirf.nirs(i).probe.(posLabelField));
                end
            end

           

            % if pos3d for field exists, check if the size matches the
            % known size, if it needs to be rotated, rotate it
            %    if nominal position field doesn't exist, use pos3D
            if(isfield(snirf.nirs(i).probe,pos3DField))
                sz=size(snirf.nirs(i).probe.(pos3DField));
                if(isnan(n))
                    sz_sel=sz(sz~=3);
                    if(~isempty(sz_sel))
                        n=sz_sel(1);
                    else
                        n=3;
                    end
                end
                
                if(any(sz==n)&&sz(1)==3&&n~=3)
                    snirf.nirs(i).probe.(pos3DField)=snirf.nirs(i).probe.(pos3DField)';
                    rotateFields=true;
                elseif(any(sz==n)&&n~=3)
                    rotateFields=false;
                elseif(~any(sz==n))
                    error('Mismatch between %s and number of elements',pos3DField);
                end

                if(~isfield(snirf.nirs(i).probe,posField))
                    snirf.nirs(i).probe.(posField)= snirf.nirs(i).probe.(pos3DField);
                end
            end

            % if pos2D for field exists, check if the size matches the
            % known size, if it needs to be rotated, rotate it
            %    if nominal position field still doesn't exist, use pos2D
            if(isfield(snirf.nirs(i).probe,pos2DField))
                sz=size(snirf.nirs(i).probe.(pos2DField));
                if(isnan(n))
                    sz_sel=sz(sz~=2);
                    if(~isempty(sz_sel))
                        n=sz_sel;
                    else
                        n=2;
                    end
                end
                
                if(any(sz==n)&&sz(1)==2&&(n~=2||rotateFields))
                    snirf.nirs(i).probe.(pos2DField)=snirf.nirs(i).probe.(pos2DField)';
                    rotateFields=true;
                elseif(any(sz==n)&&n~=2)
                    rotateFields=false;
                elseif(~any(sz==n))
                    error('Mismatch between %s and number of elements',pos2DField);
                end

                if(~isfield(snirf.nirs(i).probe,posField))
                    snirf.nirs(i).probe.(posField)= snirf.nirs(i).probe.(pos2DField);
                end

                if(rotateFields&&isfield(snirf.nirs(i).probe,pos3DField)&&n==3)
                    % Catch corner case where 2D field is rotated, but 3D
                        % field was assigned but not rotated originally
                    snirf.nirs(i).probe.(pos3DField)=snirf.nirs(i).probe.(pos3DField)';

                    if(~posFieldPreexists)
                        snirf.nirs(i).probe.(posField)=snirf.nirs(i).probe.(pos3DField);
                    end
                end
            end

            % if nominal position for field exists, check if the size matches the
            % known size, if it needs to be rotated, rotate it
            if(isfield(snirf.nirs(i).probe,posField))
                sz=size(snirf.nirs(i).probe.(posField));
                if(isnan(n))
                    sz_sel=sz(ismember(sz,[2,3]));
                    if(~isempty(sz_sel))
                        n=sz_sel(1);
                    else
                        n=3;
                        warning('Assuming %s has 3 elements',posField);
                    end
                end
                
                if((~all(sz==n)&&sz(1)~=n)||... % if n is known and mismatched rotate
                        (ismember(n,[2,3])&&(...    %if n is ambiguous 2/3
                            (~isnan(rotateFields)&&rotateFields&&~posFieldPreexists)))) % , and rotateFields is defined and true, rotate (and preexists)
                    snirf.nirs(i).probe.(posField)=snirf.nirs(i).probe.(posField)';
                elseif(ismember(n,[2,3])&&(isnan(rotateFields)||posFieldPreexists)) % if ambiguous and rotate fields not defined, or pre-existing positions
                    warning('Ambiguous probe configuration (2/3 detectors) please ensure probe coordinates are xyz as columns and element # as rows [X,Y,Z;X2,Y2,Z2]');
                elseif(~any(sz==n))
                    error('Mismatch between %s and number of elements',posField);
                end
            end

            if(~isfield(snirf.nirs(i).probe,posLabelField)&&~strcmp(posType,'landmark'))
                for j=1:n
                    snirf.nirs(i).probe.(posLabelField){j,1}=[posCapital+'-' num2str(j)];
                end
            elseif(isfield(snirf.nirs(i).probe,posLabelField)&&~strcmp(posType,'landmark'))
                snirf.nirs(i).probe.(posOrigLabelField)=snirf.nirs(i).probe.(posLabelField);
                for j=1:n
                    oldName=snirf.nirs(i).probe.(posLabelField){j};
                    newNumStr=regexprep(oldName,'\D', '');
                    if(isempty(newNumStr))
                        error('%s label must contain a number',posType);
                    end
                    snirf.nirs(i).probe.(posLabelField){j,1}=sprintf('%s-%s',posCapital,newNumStr); %NIRS toolbox requires Source-#, Detector-# format
                end
            end

            snirf.nirs(i).probe.(nField)=n;

            for j=1:n
                if(strcmp(posType,'landmark'))
                    lmName=snirf.nirs(i).probe.(posLabelField){j};
                    if(contains(lmName,'FID'))
                        Type{end+1,1}=char(snirf.nirs(i).probe.(posLabelField){j}(strfind(snirf.nirs(i).probe.(posLabelField){j},'FID'):end));
                        Name{end+1,1}=char(snirf.nirs(i).probe.(posLabelField){j}(1:strfind(snirf.nirs(i).probe.(posLabelField){j},'FID')-1));
                    else
                        Type{end+1,1}=posCapital;
                        Name{end+1,1}=lmName;
                    end 
                else
                    Type{end+1,1}=posCapital; %ex 'Source'
                    Name{end+1,1}=snirf.nirs(i).probe.(posLabelField){j};
                end
                Units{end+1,1}=snirf.nirs(i).metaDataTags.LengthUnit;
                X(end+1,1)=snirf.nirs(i).probe.(posField)(j,1);
                Y(end+1,1)=snirf.nirs(i).probe.(posField)(j,2);
                if(size(snirf.nirs(i).probe.(posField),2)==2)
                    Z(end+1,1)=0;
                else
                    Z(end+1,1)=snirf.nirs(i).probe.(posField)(j,3);
                end
            end
        end

        
        tmpdata(ii).probe.optodes_registered=table(Name,X,Y,Z,Type,Units);
        tmpdata(ii).probe.optodes=table(Name,X,Y,Z*0,Type,Units);
        
        
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

    tmpdata(ii).probe.link=table(source,detector,type);
	% add registered optodes distances to probe fixeddistances - added by Peggy Skelly
	if(isfield(snirf.nirs(i).probe,'sourcePos3D'))
		tmpdata(ii).probe.fixeddistances=tmpdata(ii).probe.swap_reg.distances;
	end
	
        for j=1:length(snirf.nirs(i).data(ii).measurementList)
            source(j,1)=snirf.nirs(i).data(ii).measurementList(j).sourceIndex;
            detector(j,1)=snirf.nirs(i).data(ii).measurementList(j).detectorIndex;
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
                
                if(~isfield(snirf.nirs(i).stim(j),'dataLabels') || isempty( snirf.nirs(i).stim(j).dataLabels))
                    snirf.nirs(i).stim(j).dataLabels={'onset','dur','amp'};
                end
                
                if(size(snirf.nirs.stim(j).data,1)==length(snirf.nirs(i).stim(j).dataLabels) & ...
                        size(snirf.nirs.stim(j).data,2)~=length(snirf.nirs(i).stim(j).dataLabels))
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

