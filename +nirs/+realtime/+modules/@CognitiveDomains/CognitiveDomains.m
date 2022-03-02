classdef CognitiveDomains < nirs.realtime.modules.AbstractModule
   properties
       model={'mindfulness.nii.gz'}
       training=30;  % number of seconds to delay the model to get an estimate of the noise 
   end
   properties(Hidden=true)
       trainingdata=0;
       projector=[];
   end
    methods
        function obj = CognitiveDomains(prevJob)
            obj.name='RT cognitive domains model';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function obj=resetThis(obj)
           obj.trainingdata=0;
           obj.projector=[];
        end
        
        function [d,t,probe,stimulus] = updateThis(obj,d,t,probe,stimulus)
            
            if(isempty(obj.projector))
            
            end
            
            
            
        end
        
    end
    
    
end