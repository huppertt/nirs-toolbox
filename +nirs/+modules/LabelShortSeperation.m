classdef LabelShortSeperation < nirs.modules.AbstractModule
    %% Label Short Seperation - adds the bool flag to the link variable
    %
    % Options:
    %     max_distance - cutoff for short distances
    
    
    properties
         max_distance = 15; % mm
    end
    
    methods
        function obj = LabelShortSeperation( prevJob )
            obj.name = 'LabelShortSeperation';
            
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            for i = 1:numel(data)
                if(ismember('ShortSeperation',data.probe.link.Properties.VariableNames))
                    data(i).probe.link.ShortSeperation=[];
                end
                ShortSeperation=(data(i).probe.distances<=obj.max_distance);
                data(i).probe.link=[data(i).probe.link table(ShortSeperation)];
            end
        end
    end
    
end

