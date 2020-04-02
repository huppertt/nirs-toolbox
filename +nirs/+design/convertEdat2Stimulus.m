function stimulus=convertEdat2Stimulus(filename)


if(isstr(filename))
    stimTable = nirs.util.readEDAT(filename);
elseif(isa(filename,'table'))
    stimTable=filename;
else
    error('unknown input type')
end

flds=stimTable.Properties.VariableNames;
isonset = false(length(flds),1);
isoffset = false(length(flds),1);
iscat = false(length(flds),1);

Names={};
for i=1:length(flds)
    if(iscellstr(stimTable.(flds{i})(1)) || isstr(stimTable.(flds{i})(1)))
        if(length(unique(stimTable.(flds{i})))>1)
            iscat(i)=true;
        end
    else
        if(~isempty(strfind(lower(flds{i}),'onsettime')))
            if(~isempty(flds{i}(1:strfind(lower(flds{i}),'_onsettime')-1)))
                isonset(i)=true;
                Names{end+1}=flds{i}(1:strfind(lower(flds{i}),'_onsettime')-1);
            end
        elseif(~isempty(strfind(lower(flds{i}),'finishtime')))
            if(~isempty(flds{i}(1:strfind(lower(flds{i}),'_finishtime')-1)))
                isoffset(i)=true;
                Names{end+1}=flds{i}(1:strfind(lower(flds{i}),'_finishtime')-1);
            end
        elseif(~isempty(strfind(lower(flds{i}),'_rttime')) | ...
                ~isempty(strfind(lower(flds{i}),'_acc')) | ...
                ~isempty(strfind(lower(flds{i}),'_rt')) | ...
                ~isempty(strfind(lower(flds{i}),'_resp')) |...
                ~isempty(strfind(lower(flds{i}),'_cresp')))
            iscat(i)=true;
        end
    end
end


Names=unique(Names);
stimulus=Dictionary;

for i=1:length(Names)
    stim=nirs.design.StimulusEvents;
    stim.name=Names{i};
    
    lst=find(isonset);
    onIdx=NaN;
    for j=1:length(lst)
        if(strcmp(lower(flds{lst(j)}),[lower(Names{i}) '_onsettime']))
            onIdx=lst(j);
            break
        end
    end
    
    
    lst=find(isoffset);
    offIdx=NaN;
    for j=1:length(lst)
        if(strcmp(lower(flds{lst(j)}),[lower(Names{i}) '_finishtime']))
            offIdx=lst(j);
            break
        end
    end
    if(~isnan(onIdx))
        for j=1:height(stimTable)
            stim.onset(j,1)=stimTable.(flds{onIdx})(j)/1000;
            if(isnan(offIdx))
                stim.dur(j,1)=1;
            else
                stim.dur(j,1)=stimTable.(flds{offIdx})(j)/1000- stim.onset(j,1);
            end
            stim.amp(j,1)=1;
            
        end
        
        stim.metadata=stimTable(:,iscat);
        
        stimulus(Names{i})=stim;
    end
    
    
    
end
