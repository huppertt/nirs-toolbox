classdef StimulusEvents
    %UNTITLED9 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        onset
        dur
        amp
    end
    
    methods 
        
        function vec = getStimVector( obj, time )
            assert( isvector( time ) )
            
            vec = zeros(size(time(:)));
            
            for i = 1:length( obj.onset )
               lst = time >= obj.onset(i) & time <= obj.onset(i)+obj.dur(i);
               vec(lst) = obj.amp(i);
            end
            
        end
        
    end
end