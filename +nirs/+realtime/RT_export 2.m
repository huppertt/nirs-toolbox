classdef RT_export < handle
    % real-time implementation of a band-pass filter
    
    properties
       datahandle=[];
    end
        
    methods

        
        function d = update(obj,d,t)
            for i=1:length(obj.datahandle)
                obj.datahandle(i).adddata(d,t);
            end
        end
            
    end
    
end