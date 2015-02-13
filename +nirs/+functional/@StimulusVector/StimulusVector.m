classdef StimulusVector
    %UNTITLED9 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        vector
        time
    end
    
    methods
        function vec = getStimVector( obj, time )
            
            vec = interp1( obj.time, obj.vector, time );

        end
    end
    
end