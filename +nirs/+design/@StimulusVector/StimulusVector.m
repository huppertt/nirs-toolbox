classdef StimulusVector
    %UNTITLED9 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        vector
        time
        regressor_no_interest=false;
        convolve_by_default=false;
    end
    
    methods
        function vec = getStimVector( obj, time )
            if(isa(obj.vector,'double'))
               vec = interp1( obj.time, obj.vector, time );
            
            else
                vec=obj.vector;
            end
            
        end
        
        function draw(obj)
            figure;
            plot(obj.time,obj.vector);
            xlabel('Time (sec)')
            ylabel(obj.name);
        end
        
    end
    
end