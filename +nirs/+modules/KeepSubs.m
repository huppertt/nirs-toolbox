classdef KeepSubs < nirs.modules.AbstractModule
    %% KeepSubs - Removes all subjects except those specified.
    %
    % Options:
    %     listOfFilters - [N x 2] cell array of demographics key-value pairs
    %                     to keep (e.g. {'Group', 'Control'; 'Gender', 'M' }
    properties
        listOfFilters = {}; % cell array of demographic criteria to match
    end
    
    methods
        function obj = KeepSubs( prevJob )
            if nargin > 0; obj.prevJob = prevJob; end
            obj.name = 'Keep Just These Subjects';
        end
        
        function data = runThis( obj, data )
            
            if isempty(obj.listOfFilters), warning('No subject filter specified. Skipping.'); return; end
            
            keepers = true(size(data));
            
            for i = 1:length(data)
                
                demo = data(i).demographics;
                
                for j = 1:size(obj.listOfFilters,1)
                
                    if ~strcmpi( demo( obj.listOfFilters{j,1} ) , obj.listOfFilters{j,2} )
                        
                        keepers(i) = false;
                        
                    end
                end
            end
            
            data = data(keepers);
        end
    end
end

