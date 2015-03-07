classdef RenameStims < nirs.functional.AbstractModule
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        list = {};
    end
    
    methods
        function obj = RenameStims( prevJob )
           obj.name = 'Rename Stim Conditions';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = execute( obj, data )

            [names, idx] = nirs.functional.getStimNames( data );
            for i = 1:length(obj.list)
                lst = strcmp(names, obj.list{i,1});
                names(lst) = repmat( (obj.list(i,2)), [sum(lst) 1] );
            end
            
            for i = 1:length(data)
                lst = idx == i;
                
                if any(lst)
                    keys = names(lst);
                    values = data(i).stimulus.values;

                    ukeys = unique(keys,'stable');
                    for j = 1:length(ukeys)
                       lst = strcmp(ukeys{j},keys);
                       uvalues{j} = nirs.functional.mergeStims( values(lst), ukeys{j} );
                    end

                    data(i).stimulus = nirs.HashTable(ukeys, uvalues);
                end
            end
            
            
        end
        

        
        function options = getOptions( obj )
            options = [];
        end
           
        function obj = putOptions( obj, options )
        end
    end
    
end

