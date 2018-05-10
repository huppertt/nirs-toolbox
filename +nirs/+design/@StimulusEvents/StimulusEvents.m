classdef StimulusEvents
    properties
        name
        onset
        dur
        amp
        metadata=table;
    end
    
    methods
        function obj = StimulusEvents( name, onset, dur, amp )
           if nargin > 0, obj.name  = name; end
           if nargin > 1, obj.onset = onset; end
           if nargin > 2, obj.dur   = dur; end
           if nargin > 3, obj.amp   = amp; end
           
        end
        
%         function disp(obj)
%             disp(['nirs.design.StimulusEvents [' num2str(length(obj.onset)) ' events']);
%         end
%         
        
        function vec = getStimVector( obj, time )
            assert( isvector( time ) )
            
            vec = zeros(size(time(:)));
            
            for i = 1:length( obj.onset )
                % list of points for this onset
                lst = time >= obj.onset(i) & time < obj.onset(i)+obj.dur(i);
                
                if(isempty(find(lst)))
                    % for really short durations this can happen because of
                    % the >= vs < issue
                    lst=min(find(time >= obj.onset(i)));
                end
                
                % set them to correct amplitude
                vec(lst) = vec(lst) + obj.amp(i);
                    
            end
            
        end
        
    end
end