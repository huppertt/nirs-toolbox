classdef DiscardStimWrongCount < nirs.modules.AbstractModule
    %% DiscardStimWrongCount - Removes conditions w/ incorrect number of stimuli.
    %
    % Options:
    %     listOfCounts - [N x 2] cell array of condition name and count pairs
    %                     to keep (e.g. {'GoNoGo', 16; 'WM', 4 }
    properties
        listOfCounts = {}; % cell array of stimulus counts for each condition
    end
    
    methods
        function obj = DiscardStimWrongCount( prevJob )
            if nargin > 0; obj.prevJob = prevJob; end
            obj.name = 'Remove conditions with wrong number of stimuli';
        end
        
        function data = runThis( obj, data )
            
            if isempty(obj.listOfCounts), warning('No list of stimulus counts specified. Skipping.'); return; end
                       
            for i = 1:length(data)
                
                stimnames = data(i).stimulus.keys;
                
                for j = 1:size(obj.listOfCounts,1)
                
                    if any( strcmpi( stimnames , obj.listOfCounts{j,1} ) )
                    
                        stim = data(i).stimulus( obj.listOfCounts{j,1} );
                        
                        if length(stim.onset) ~= obj.listOfCounts{j,2}
                            
                            data(i).stimulus = data(i).stimulus.remove( obj.listOfCounts{j,1} );
                            
                        end
                    end
                end
            end          
        end
    end
end

