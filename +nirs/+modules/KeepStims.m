classdef KeepStims < nirs.modules.AbstractModule
    %% KeepStims - Removes all stims except those specified.
    %
    % Options:
    %     listOfStims - cell array of stim names
    properties
        listOfStims = {}; % cell array of stim names
    end
    
    methods
        function obj = KeepStims( prevJob )
            if nargin > 0; obj.prevJob = prevJob; end
            obj.name = 'Keep Just These Stims';
        end
        
        function data = runThis( obj, data )
            for i = 1:length(data)
                if(isa(data(i),'nirs.core.Data') | isa(data(i),'eeg.core.Data') | isa(data(i),'dtseries.core.Data'))
                    stim = data(i).stimulus;
                    
                    % get values for stims we keep
                    keys = obj.listOfStims;
                    vals = stim( keys );
                    
                    if ~iscell(vals)
                        vals = {vals};
                    end
                    
                    % lst of nonempty
                    lst = ~cellfun(@isempty, vals);
                    
                    % new stimulus
                    data(i).stimulus = Dictionary(keys(lst), vals(lst));
                    
                elseif(isa(data(i),'nirs.core.sFCStats'))
                    cond=data(i).conditions;
                    lst=find(ismember(cond,obj.listOfStims));
                    data(i).R=data(i).R(:,:,lst);
                    data(i).dfe=data(i).dfe(lst);
                    data(i).conditions={data(i).conditions{lst}};
                else
                    cond=data(i).conditions;
                    cond={cond{ismember(cond,obj.listOfStims)}};
                    data(i)=data(i).ttest(cond);
                end
            end
            
            if isempty(nirs.getStimNames(data))
                warning('No stim conditions left. Did you provide the correct stim names?')
                
            end
        end
    end
end

