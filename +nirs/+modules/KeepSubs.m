classdef KeepSubs < nirs.modules.AbstractModule
    %% KeepSubs - Removes all subjects except those specified.
    %
    % Options:
    %     listOfFilters - [N x 2] cell array of demographics key-value pairs
    %                     to keep (e.g. {'Group', 'Control'; 'Gender', 'M'; 'SubjectID', {'Sub01','Sub03','Sub09'} }
    %     remove_matches - logical (true inverts selection criteria; default: false)
    properties
        listOfFilters = {}; % cell array of demographic criteria to match
        remove_matches = false; % Inverts the selection criteria to remove those who match
    end
    
    methods
        function obj = KeepSubs( prevJob )
            if nargin > 0; obj.prevJob = prevJob; end
            obj.name = 'Keep Just These Subjects';
        end
        
        function data = runThis( obj, data )
            
            if isempty(obj.listOfFilters), warning('No subject filter specified. Skipping.'); return; end
            
            keepers = true(size(data));
            
            for i = 1:numel(data)
                
                demo = data(i).demographics;
                
                for j = 1:size(obj.listOfFilters,1)
                
                    value = demo( obj.listOfFilters{j,1} );
                    allowable = obj.listOfFilters{j,2};
                    is_match = any(strcmpi( allowable , value  ));
                    
                    if ~obj.remove_matches && ~is_match
                        keepers(i) = false;
                    end
                    if obj.remove_matches && is_match
                        keepers(i) = false;
                    end
                    
                end
            end
            
            data = data(keepers);
        end
    end
end

