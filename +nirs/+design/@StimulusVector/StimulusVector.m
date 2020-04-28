classdef StimulusVector
    %UNTITLED9 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        vector
        time
<<<<<<< HEAD
=======
        regressor_no_interest=false;
>>>>>>> e5f3a84a411a47d49ab3f360371dbfa4180f0a9f
    end
    
    methods
        function vec = getStimVector( obj, time )
<<<<<<< HEAD
            
            vec = interp1( obj.time, obj.vector, time );

=======
        vec = interp1( obj.time, obj.vector, time );
>>>>>>> e5f3a84a411a47d49ab3f360371dbfa4180f0a9f
        end
        
        function draw(obj)
            figure;
            plot(obj.time,obj.vector);
            xlabel('Time (sec)')
            ylabel(obj.name);
        end
        
    end
    
end