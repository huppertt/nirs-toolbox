classdef DiscardStims < nirs.modules.AbstractModule
%% DiscardStims - Removes specified stim conditions from design.
% 
% Options: 
%     listOfStims - cell array of stim names

    properties
        listOfStims = {}; % cell array of stim names
    end
    
    methods
        function obj = DiscardStims( prevJob )
         	if nargin > 0, obj.prevJob = prevJob; end
            obj.name = 'Discard Stim Conditions';
        end
        
        function data = runThis( obj, data )
            % get current stim names
            stimNames = nirs.getStimNames(data);

            % loop through and remove stims
            for i = 1:length(data)
                data(i).stimulus = data(i).stimulus.remove( obj.listOfStims );
            end

            % compare old and new list of stims
            if length(stimNames) == length(nirs.getStimNames(data))
                warning( 'No stims discarded.  Did you provide the correct stim names?' )
            end
        end
    end
    
end

