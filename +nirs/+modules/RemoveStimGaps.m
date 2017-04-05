classdef RemoveStimGaps < nirs.modules.AbstractModule
%% RemoveStimGaps - Merges task blocks that are only separated by a user-specified interval.
% 
% Options: 
%     maxDuration - scalar maximum duration (in samples) of gap to remove

    properties
        maxDuration = 1; % maximum duration of gap to remove
    end
    
    methods
        function obj = RemoveStimGaps( prevJob )
         	if nargin > 0, obj.prevJob = prevJob; end
            obj.name = 'Remove gaps in stim block';
        end
        
        function data = runThis( obj, data )
            
            assert(isa(data,'nirs.core.Data'));
            
            for i = 1:length(data)
                stimNames = nirs.getStimNames(data(i));
                for j = 1:length(stimNames)
                    
                    stims = data(i).stimulus(stimNames{j});
                    onset = stims.onset * data(i).Fs;
                    dur = stims.dur * data(i).Fs;
                    amp = stims.amp;
                    offset = onset + dur;
                    diff = onset(2:end) - offset(1:end-1);
                    mergeInds = find(diff<=obj.maxDuration);
                    
                    if ~isempty(mergeInds)
                    
                        for k = length(mergeInds):-1:1

                            onset(mergeInds(k)+1) = [];
                            offset(mergeInds(k)) = [];
                            amp(mergeInds(k)+1) = [];

                        end
                        stims.onset = onset ./ data(i).Fs;
                        stims.dur = (offset - onset) ./ data(i).Fs;
                        stims.amp = amp;

                        data(i).stimulus(stimNames{j}) = stims;
                        
                    end
                end
                
            end            
        end
    end
    
end

