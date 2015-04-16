classdef DiscardStims < nirs.modules.AbstractModule
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        listOfStims     = {};
    end
    
    properties ( Hidden = true )
        keepOrDiscard   = 'discard';
    end
    
    methods
        function obj = DiscardStims( prevJob )
           obj.name = 'Discard Stim Conditions';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            
            names = unique( nirs.getStimNames(data) );
            
            % find all in names that match the list of stims
            lst = zeros(size(names));
            for i = 1:length(obj.listOfStims)
                lst = lst  | strcmp(obj.listOfStims{i}, names);
            end
            
            %if flag == discard invert list
            if strcmp(obj.keepOrDiscard,'keep')
                lst = ~lst;
            end
            
            names = names(lst);
            
            if isempty(names)
                warning( 'No stims discarded.  Did you provide the correct stim names?' )
            else
                % loop over files and discard stims
                for i = 1:length(data)
                   data(i).stimulus = data(i).stimulus.delete( names );
                end
            end

%             % get all stim names and data file indices
%             [names, idx] = nirs.getStimNames( data );
%             
%             % find all in names that match obj.list
%             lst = zeros(size(idx));
%             for i = 1:length(obj.listOfStims)
%                 lst = lst  | strcmp(names, obj.listOfStims{i});
%             end
%             
%             % if flag == discard invert list
%             if strcmp(obj.keepOrDiscard,'discard')
%                 lst = ~lst;
%             end
%             
%             % loop through files and only keep stims in lst
%             for i = 1:length(data)
%                 keys    = data(i).stimulus.keys;
%                 values  = data(i).stimulus.values;
%                 
%                 keys    = keys( lst(i==idx) );
%                 values  = values( lst(i==idx) );
%                 
%                 data(i).stimulus = Dictionary(keys, values);
%             end
        end
    end
    
end

