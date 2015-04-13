classdef RemoveStimless < nirs.modules.AbstractModule

    methods
        function obj = RemoveStimless( prevJob )
           obj.name = 'Remove Files w/o Stim';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:length(data)
                if isempty(data(i).stimulus.keys)
                    lst(i) = false;
                else
                    lst(i) = true;
                end
            end
            
            data = data(lst);
            
        end
    end
    
end

