classdef DiscardStimWrongCount < nirs.modules.AbstractModule
    %% DiscardStimWrongCount - Removes conditions w/ incorrect number of stimuli.
    %
    % Options:
    %     listOfCounts - [N x 2] cell array of condition name and count pairs
    %                     to keep (e.g. {'GoNoGo', 16; 'WM', 4 }
    %                       OR
    %                    [N x 3] cell array including comparison operator
    %                      (e.g. {'GoNoGo', '<=', 16;   'GoNoGo', '>', 8 };
    %                      would remove GoNoGo if # was less than 9 or
    %                      greater than 16
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
            
            if size(obj.listOfCounts,2)<3
                obj.listOfCounts(:,3) = obj.listOfCounts(:,2);
                obj.listOfCounts(:,2) = {'=='};
            end
            
            for i = 1:numel(data)
                
                stimnames = data(i).stimulus.keys;
                
                for j = 1:size(obj.listOfCounts,1)
                
                    if any( strcmpi( stimnames , obj.listOfCounts{j,1} ) )
                    
                        stim = data(i).stimulus( obj.listOfCounts{j,1} );
                        if isempty(stim), continue; end
                        
                        test = sprintf('length(stim.onset) %s obj.listOfCounts{j,3}',obj.listOfCounts{j,2});
                        
                        if ~eval(test)
                            
                            data(i).stimulus = data(i).stimulus.remove( obj.listOfCounts{j,1} );
                            
                        end
                    end
                end
            end          
        end
    end
end

