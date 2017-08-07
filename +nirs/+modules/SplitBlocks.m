classdef SplitBlocks < nirs.modules.AbstractModule
%% SplitBlocks - Splits long task blocks into multiple smaller blocks
% 
    properties
        maxtime=30;  % time in seconds
    end
    methods
        function obj = SplitBlocks( prevJob )
           obj.name = 'Remove Short Stims';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:length(data)
                keys = data(i).stimulus.keys;
                for j = 1:length(keys)
                    stim = data(i).stimulus(keys{j});
                    stim2 = nirs.design.StimulusEvents(keys{j});
                    for k = 1:length(stim.onset) % Loop over actual task blocks
                        numblocks = ceil(stim.dur(k)/obj.maxtime);
                        stim2.onset(end+1) = stim.onset(k);
                        stim2.dur(end+1) = min(stim.dur(k),obj.maxtime);
                        stim2.amp(end+1) = stim.amp(k);
                        offset = stim.onset(k) + stim.dur(k);
                        for l = 2:numblocks % Loop over each new block segment
                            stim2.onset(end+1) = stim2.onset(end)+stim2.dur(end);
                            stim2.dur(end+1) = min(offset-stim2.onset(end),obj.maxtime);
                            stim2.amp(end+1) = stim.amp(k);
                        end
                    end
                    data(i).stimulus(keys{j}) = stim2;                    
                end
            end
        end
    end
    
end

