classdef LabelnotEnoughDistance < nirs.modules.AbstractModule
    %% Label too Long Distance - adds the bool flag to the link variable
    %
    % Options:
    %     min_distance - cutoff for min distance
    %     max_distance - cutoff for max distance
    
    
    properties
         min_distance = 10; % mm
         max_distance = 25; % mm
    end
    
    methods
        function obj = LabeltooLongSeperation( prevJob )
            obj.name = 'LabelnotEnoughDistance';
            
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            for i = 1:numel(data)
                if(ismember('notEnoughDistance',data(i).probe.link.Properties.VariableNames))
                    data(i).probe.link.notEnoughDistance=[];
                end
                notEnoughDistance=and(data(i).probe.distances>=obj.min_distance, ...
                    data(i).probe.distances<=obj.max_distance);
                data(i).probe.link=[data(i).probe.link table(notEnoughDistance)];
            end
        end
    end
    
end

