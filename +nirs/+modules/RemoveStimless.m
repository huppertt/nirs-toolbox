classdef RemoveStimless < nirs.modules.AbstractModule
%% RemoveStimless - Removes files with no stimulus information.
% 


    methods
        function obj = RemoveStimless( prevJob )
           obj.name = 'Remove Files w/o Stim';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:length(data)
                if(isa(data(i),'nirs.core.Data') | isa(data(i),'eeg.core.Data'))
                    if isempty(data(i).stimulus.keys)
                        lst(i) = false;
                    else
                        lst(i) = true;
                    end
                else
                    if isempty(data(i).conditions)
                        lst(i) = false;
                    else
                        lst(i) = true;
                    end
                end
            end
            
            data = data(lst);
            
        end
    end
    
end

