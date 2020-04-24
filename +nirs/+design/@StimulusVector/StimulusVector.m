classdef StimulusVector
    %UNTITLED9 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        vector
        time
        regressor_no_interest=false;
    end
    
    methods
        function vec = getStimVector( obj, time )
        vec = interp1( obj.time, obj.vector, time );
        end
        
        function draw(obj)
            figure;
            plot(obj.time,obj.vector);
            xlabel('Time (sec)')
            ylabel(obj.name);
        end
        
    end
    
end