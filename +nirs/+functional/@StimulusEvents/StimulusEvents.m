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
        function obj = StimulusEvents( name, onset, dur, amp )
           if nargin > 0, obj.name = name; end
           if nargin > 1, obj.onset = onset; end
           if nargin > 2, obj.dur = dur; end
           if nargin > 3, obj.amp = amp; end
        end
        
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