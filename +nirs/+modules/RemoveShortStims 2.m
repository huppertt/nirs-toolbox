classdef RemoveShortStims< nirs.modules.AbstractModule
%% RemoveShortStims - Removes stimuli with shorter than minimum time in length
% 
    properties
        mintime=30;  % time in seconds
    end
    methods
        function obj = RemoveShortStims( prevJob )
           obj.name = 'Remove Short Stims';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:numel(data)
                keys = data(i).stimulus.keys;
                for j = 1:length(keys)
                    stim = data(i).stimulus(keys{j});
                    is_short = (stim.dur < obj.mintime);
                    if all(is_short)
                        data(i).stimulus = data(i).stimulus.remove(keys{j});
                    elseif any(is_short)
                        stim.onset(is_short) = [];
                        stim.dur(is_short) = [];
                        stim.amp(is_short) = [];
                        data(i).stimulus(keys{j}) = stim;
                    end
                end
            end
        end
    end
    
end

