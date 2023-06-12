classdef RenameStims < nirs.modules.AbstractModule
%% RenameStims - Renames stimulus conditions.
% 
% Options:
%     listOfChanges - % n x 2 cell array of changes
%     
% Example:
%     j = nirs.modules.RenameStims();
%     j.listOfChanges = { ...
%         'miSSpellledHOrror',        'NiceName1';
%         'HorRRRRORTWO',             'NiceName2';
%         'OMGThisNameIsTooDangLong', 'NiceName3'
%         };
%
%     j.run( raw );
%     
% Notes:
%     If you rename a condition to one that exists, the two stims are merged.
    
    properties
        listOfChanges = {}; % n x 2 cell array of changes
    end
    
    methods
        function obj = RenameStims( prevJob )
           obj.name = 'Rename Stim Conditions';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            
            % get all stim names across all files
            [names, idx] = nirs.getStimNames( data );
            for i = 1:size(obj.listOfChanges,1)
                lst = strcmp(names, obj.listOfChanges{i,1});
                names(lst) = repmat( (obj.listOfChanges(i,2)), [sum(lst) 1] );
            end
            
            % for each file rename stims
            for i = 1:numel(data)
                lst = idx == i;
                
                if any(lst)
                    if(isa(data(i),'nirs.core.Data') || isa(data(i),'eeg.core.Data') ||...
                            isa(data(i),'nirs.core.GenericData'))
                    keys = names(lst)';
                    values = data(i).stimulus.values;

                    ukeys = unique(keys,'stable');
                    uvals = {};
                    for j = 1:length(ukeys)
                       lst = strcmp(ukeys{j},keys);
                       uvals{j} = nirs.design.mergeStims( values(lst), ukeys{j} );
                    end

                    data(i).stimulus = Dictionary(ukeys, uvals);
                    elseif(isa(data(i),'nirs.core.ChannelStats') || isa(data(i),'eeg.core.ChannelStats') )
                        [loca,locb]=ismember(data(i).variables.cond,{obj.listOfChanges{:,1}}); 
                        for jj=1:length(loca)
                            if(loca(jj))
                                data(i).variables.cond{jj}=obj.listOfChanges{locb(jj),2};
                            end
                        end
                        
                        try
                            keys=data(i).basis.stim.keys;
                            s=Dictionary;
                            for jj=1:length(keys)
                                lst=find(ismember({obj.listOfChanges{:,1}},keys{jj}));
                                if(isempty(lst))
                                    s(keys{jj})=data(i).basis.stim(keys{jj});
                                else
                                    s(obj.listOfChanges{lst,2})=data(i).basis.stim(keys{jj});
                                end
                            end
                            data(i).basis.stim=s;
                        end
                    end
                end
            end
        end
    end
    
end

