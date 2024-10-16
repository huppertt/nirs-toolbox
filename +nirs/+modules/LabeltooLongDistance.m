classdef LabeltooLongDistance < nirs.modules.AbstractModule
    %% Label too Long Distance - adds the bool flag to the link variable
    %
    % Options:
    %     min_distance - cutoff for too long distance
    
    
    properties
         min_distance = 35; % mm
    end
    
    methods
        function obj = LabeltooLongSeperation( prevJob )
            obj.name = 'LabeltooLongSeperation';
            
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            for i = 1:numel(data)
                if(ismember('tooLongDistance',data(i).probe.link.Properties.VariableNames))
                    data(i).probe.link.tooLongDistance=[];
                end
                tooLongDistance=(data(i).probe.distances>=obj.min_distance);
                data(i).probe.link=[data(i).probe.link table(tooLongDistance)];
            end
        end
    end
    
end

