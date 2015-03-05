classdef DiscardStims < nirs.functional.AbstractModule
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        list = {};
        flag = 'discard';
    end
    
    methods
        function obj = RenameStims( prevJob )
           obj.name = 'Discard Stim Conditions';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = execute( obj, data )

            % get all stim names and data file indices
            [names, idx] = nirs.functional.getStimNames( data );
            
            % find all in names that match obj.list
            lst = zeros(size(idx));
            for i = 1:length(obj.list)
                lst = lst  | strcmp(names, obj.list{i});
            end
            
            % if flag == discard invert list
            if strcmp(obj.flag,'discard')
                lst = ~lst;
            end
            
            % loop through files and only keep stims in lst
            for i = 1:length(data)
                keys    = data(i).stimulus.keys;
                values  = data(i).stimulus.values;
                
                keys    = keys( lst(i==idx) );
                values  = values( lst(i==idx) );
                
                data(i).stimulus = nirs.HashTable(keys, values);
            end
        end
        
        function options = getOptions( obj )
            options = [];
        end
           
        function obj = putOptions( obj, options )
        end
    end
    
end

