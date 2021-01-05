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
                          
        regex  = false;   % logical dictating whether to match condition names based on a regular experession
                          % Ex. '^[012]back$' would match 0back, 1back, 2back
                          % Ex. ':01$' would keep anything ending in :01
                          % default: false
    end
    
    methods
        function obj = KeepStims( prevJob )
            if nargin > 0; obj.prevJob = prevJob; end
            obj.name = 'Keep Just These Stims';
        end
        
        function data = runThis( obj, data )
            
            if obj.regex
                checkmatch = @(conds,reference) ~cellfun(@isempty,nirs.util.regexpi(conds,reference,'forcecelloutput'));
            else
                checkmatch = @(conds,reference) strcmpi(conds,reference);
            end
            
            obj.required = logical(obj.required);
            if numel(obj.required)==1
                obj.required = repmat(obj.required,size(obj.listOfStims));
            end
            if length(obj.required) ~= length(obj.listOfStims)
                error('List of stims and required paramters not compatible. Lengths: %i and %i',length(obj.required),length(obj.listOfStims));
            end
            
            bad_subjects = false(size(data));
            
            for i = 1:numel(data)
                
                if(isa(data(i),'nirs.core.Data') || isa(data(i),'eeg.core.Data') || isa(data(i),'dtseries.core.Data'))
                    conds=data(i).stimulus.keys;
                    isgoodcond=false(size(conds));
                    for j=1:length(obj.listOfStims)
                        match = checkmatch(conds,obj.listOfStims{j});
                        isgoodcond = isgoodcond | match;
                        if obj.required(j) && ~any(match)
                            bad_subjects(i)=true;
                        end
                    end
                    
                    badconds = conds(~isgoodcond);
                    for j=1:length(badconds)
                        data(i).stimulus = data(i).stimulus.remove(badconds{j});
                    end
                    
                elseif(isa(data(i),'nirs.core.sFCStats'))
                    conds=data(i).conditions;
                    isgoodcond=false(size(conds));
                    for j=1:length(obj.listOfStims)
                        match = checkmatch(conds,obj.listOfStims{j});
                        isgoodcond = isgoodcond | match;
                        if obj.required(j) && ~any(match)
                            bad_subjects(i)=true;
                        end
                    end
                    data(i).R=data(i).R(:,:,isgoodcond);
                    data(i).dfe=data(i).dfe(isgoodcond);
                    data(i).conditions={data(i).conditions{isgoodcond}};
                    
                elseif(isa(data(i),'nirs.core.ChannelStats'))
                    conds=data(i).variables.cond;
                    isgoodcond=false(size(conds));
                    for j=1:length(obj.listOfStims)
                        match = checkmatch(conds,obj.listOfStims{j});
                        isgoodcond = isgoodcond | match;
                        if obj.required(j) && ~match
                            bad_subjects(i)=true;
                        end
                    end
                    data(i).variables(~isgoodcond,:)=[];
                    data(i).beta(~isgoodcond) = [];
                    data(i).covb(~isgoodcond,:) = [];
                    data(i).covb(:,~isgoodcond) = [];
                    
                elseif(isa(data(i),'nirs.core.ChannelFStats'))
                    conds=data(i).variables.cond;
                    isgoodcond=false(size(conds));
                    for j=1:length(obj.listOfStims)
                        match = checkmatch(conds,obj.listOfStims{j});
                        isgoodcond = isgoodcond | match;
                        if obj.required(j) && ~match
                            bad_subjects(i)=true;
                        end
                    end
                    data(i).variables(~isgoodcond,:)=[];
                    data(i).F(~isgoodcond) = [];
                    data(i).df1(~isgoodcond) = [];
                    data(i).df2(~isgoodcond) = [];
                    if any(~ismember(obj.listOfStims(obj.required),conds(isgoodcond)))
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

