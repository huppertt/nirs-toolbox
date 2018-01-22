classdef KeepStims < nirs.modules.AbstractModule
    %% KeepStims - Removes all stims except those specified.
    %
    % Options:
    %     listOfStims - cell array of stim names
    properties
        listOfStims = {}; % cell array of stim names
        required = false; % logical scalar or array [1 x ncond] indicating
                          % whether stimuli are required (subjects missing them will be removed)
                          % default: false
    end
    
    methods
        function obj = KeepStims( prevJob )
            if nargin > 0; obj.prevJob = prevJob; end
            obj.name = 'Keep Just These Stims';
        end
        
        function data = runThis( obj, data )
            
            obj.required = logical(obj.required);
            if numel(obj.required)==1
                obj.required = repmat(obj.required,size(obj.listOfStims));
            end
            if length(obj.required) ~= length(obj.listOfStims)
                error('List of stims and required paramters not compatible. Lengths: %i and %i',length(obj.required),length(obj.listOfStims));
            end
            bad_subjects = false(size(data));
            
            for i = 1:numel(data)
                
                if(isa(data(i),'nirs.core.Data') | isa(data(i),'eeg.core.Data') | isa(data(i),'dtseries.core.Data'))
                    stim = data(i).stimulus;
                    
                    % get values for stims we keep
                    keys = obj.listOfStims;
                    if(isempty(keys))
                        data(i).stimulus=Dictionary;
                        continue;
                    end
                    
                    vals = stim( keys );
                    
                    if ~iscell(vals)
                        vals = {vals};
                    end
                    
                    % lst of nonempty
                    lst = ~cellfun(@isempty, vals);
                    
                    % new stimulus
                    data(i).stimulus = Dictionary(keys(lst), vals(lst));
                    
                    if any( obj.required(:) & ~lst(:) )
                        bad_subjects(i) = true;
                    end
                    
                elseif(isa(data(i),'nirs.core.sFCStats'))
                    cond=data(i).conditions;
                    lst=find(ismember(cond,obj.listOfStims));
                    data(i).R=data(i).R(:,:,lst);
                    data(i).dfe=data(i).dfe(lst);
                    data(i).conditions={data(i).conditions{lst}};
                    if any(~ismember(obj.listOfStims(obj.required),cond))
                        bad_subjects(i) = 1;
                    end
                    
                elseif(isa(data(i),'nirs.core.ChannelStats'))
                    conds=data(i).variables.cond;
                    isgood=false(size(conds));
                    for j=1:length(obj.listOfStims)
                        isgood = isgood | strcmp(conds,obj.listOfStims{j});
                    end
                    data(i).variables(~isgood,:)=[];
                    data(i).beta(~isgood) = [];
                    data(i).covb(~isgood,:) = [];
                    data(i).covb(:,~isgood) = [];
                    if any(~ismember(obj.listOfStims(obj.required),conds(isgood)))
                        bad_subjects(i) = 1;
                    end
                    
                elseif(isa(data(i),'nirs.core.ChannelFStats'))
                    conds=data(i).variables.cond;
                    isgood=false(size(conds));
                    for j=1:length(obj.listOfStims)
                        isgood = isgood | strcmp(conds,obj.listOfStims{j});
                    end
                    data(i).variables(~isgood,:)=[];
                    data(i).F(~isgood) = [];
                    data(i).df1(~isgood) = [];
                    data(i).df2(~isgood) = [];
                    if any(~ismember(obj.listOfStims(obj.required),conds(isgood)))
                        bad_subjects(i) = 1;
                    end
                    
                else
                    cond=data(i).conditions;
                    cond={cond{ismember(cond,obj.listOfStims)}};
                    data(i)=data(i).ttest(cond);
                    if any(~ismember(obj.listOfStims(obj.required),cond))
                        bad_subjects(i) = 1;
                    end
                end
            end
            
            if isempty(nirs.getStimNames(data))
                warning('No stim conditions left. Did you provide the correct stim names?')
                
            end
            
            data(bad_subjects) = [];
        end
    end
end

