classdef RemoveStimGaps < nirs.modules.AbstractModule
    %% RemoveStimGaps - Merges task blocks that are only separated by a user-specified interval.
    %
    % Options:
    %     maxDuration - scalar maximum duration (in samples) of gap to remove
    %     minISI - minumum inter-stumulus interval
    % Note- one or both methods can be used.  Set the unused option to NaN;
    
    properties
        maxDuration = 1; % maximum duration of gap to remove
        minISI =NaN;
    end
    
    methods
        function obj = RemoveStimGaps( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            obj.name = 'Remove gaps in stim block';
        end
        
        function data = runThis( obj, data )
            
            if(isnan(obj.minISI))
                obj.minISI=0;
            end
            if(isnan(obj.maxDuration))
                obj.maxDuration=Inf;
            end
            
            
            assert(isa(data,'nirs.core.Data'));
            
            for i = 1:length(data)
                
                merge2previous=[];
                mergeFrom=[];
                stimNames = nirs.getStimNames(data(i));
                
                onsets=[]; offsets=[];  stimtypes=[];
                for j = 1:length(stimNames)
                    stims = data(i).stimulus(stimNames{j});
                    onsets=[onsets; stims.onset(:)];
                    offsets=[offsets; stims.onset(:)+stims.dur(:)];
                    stimtypes=[stimtypes; ones(size(stims.onset(:)))*j];
                end
                [onsets, sIdx]=sort(onsets);
                offsets=offsets(sIdx);
                stimtypes=stimtypes(sIdx);
                
                for j = 1:length(stimNames)
                    stims = data(i).stimulus(stimNames{j});
                    for k=1:length(stims.onset)
                        os=stims.onset(k);
                        thisIdx=find(onsets==os & stimtypes==j);
                        
                        IdxLastStim=max(find(onsets<os));
                        timefromlaststim = onsets(IdxLastStim);
                        ISI2fromlaststim = os-timefromlaststim;
                        gaptimefromlaststim = os-offsets(IdxLastStim);
                        gapsampfromlaststim = gaptimefromlaststim*data(i).Fs;
                        
                        merge2previous(thisIdx)=false;
                        mergewith(thisIdx)=NaN;
                        if(~isempty(IdxLastStim))
                            if(gapsampfromlaststim<=obj.maxDuration || ISI2fromlaststim<=obj.minISI)
                                merge2previous(thisIdx)=true;
                                mergewith(thisIdx)=IdxLastStim;
                            end
                            
                            %else- this is the first event in the file
                        end
                        
                    end
                end
                
                lst=find(merge2previous);
                
                for j=1:length(lst)
                    mergeFrom(j,1) = stimtypes(lst(j));
                    st=data(i).stimulus(stimNames{mergeFrom(j,1)});
                    mergeFrom(j,2) = find(st.onset==onsets(lst(j)));
                    
                    mergeTo(j,1) = stimtypes(mergewith(lst(j)));
                    st=data(i).stimulus(stimNames{mergeTo(j,1)});
                    mergeTo(j,2) = find(st.onset==onsets(mergewith(lst(j))));
                end
                
                for j=1:length(lst)
                    stTo=data(i).stimulus(stimNames{mergeTo(j,1)});
                    stFrom=data(i).stimulus(stimNames{mergeFrom(j,1)});
                    stTo.dur(mergeTo(j,2))=stFrom.onset(mergeFrom(j,2))+stFrom.dur(mergeFrom(j,2))-stTo.onset(mergeTo(j,2));
                  
                    data(i).stimulus(stimNames{mergeTo(j,1)})=stTo;
                end
                if(~isempty(lst))
                    for j=1:length(stimNames)
                        lst2=mergeFrom(find(mergeFrom(:,1)==j),2);
                        
                        stFrom=data(i).stimulus(stimNames{j});
                        stFrom.onset(lst2)=[];
                        stFrom.dur(lst2)=[];
                        stFrom.amp(lst2)=[];
                        data(i).stimulus(stimNames{j})=stFrom;
                    end
                end
                
            end
        end
    end
end
