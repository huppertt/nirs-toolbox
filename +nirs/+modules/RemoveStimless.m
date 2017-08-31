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
            lst = false(size(data));
            for i = 1:length(data)
                if(isa(data(i),'nirs.core.Data') | isa(data(i),'eeg.core.Data'))
                    if ~isempty(data(i).stimulus.keys)
                        lst(i) = true;
                        
                        on=[];
                        for j=1:data(i).stimulus.count
                            ss=data(i).stimulus.values{j};
                            on=[on; ss.onset(:)];
                        end
                        if(all(on> max(data(i).time)))
                            lst(i)=false;
                        end
                    end
                else
                    if ~isempty(data(i).conditions)
                        lst(i) = true;
                    end
                end
            end
            
            data = data(lst);
            
        end
    end
    
end

