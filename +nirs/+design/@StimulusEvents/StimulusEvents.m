classdef StimulusEvents
    properties
        name
        onset
        dur
        amp
    end
    
    methods
        function obj = StimulusEvents( name, onset, dur, amp )
           if nargin > 0, obj.name  = name; end
           if nargin > 1, obj.onset = onset; end
           if nargin > 2, obj.dur   = dur; end
           if nargin > 3, obj.amp   = amp; end
        end
        
        function vec = getStimVector( obj, time )
            assert( isvector( time ) )
            
            vec = zeros(size(time(:)));
            
            for i = 1:length( obj.onset )
                % list of points for this onset
                lst = time >= obj.onset(i) & time <= obj.onset(i)+obj.dur(i);
                
                % set them to correct amplitude
                vec(lst) = vec(lst) + obj.amp(i);
            end
            
        end
        
    end
end