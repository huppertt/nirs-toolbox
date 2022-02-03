classdef TrimBaseline < nirs.modules.AbstractModule
%% TrimBaseline - Removes excessive baseline at the beginning or end of file.
%
% Options:
%     preBaseline  - maximum baseline (seconds) at the beginning of scan
%     postBaseline - maximum baseline (seconds) after final task period ends
%     resetTime    - flag indicating whether the time vector should be reset so t(1)=0 (default: false)
    
    properties
        preBaseline  = 30;  % maximum baseline (seconds) at the beginning of scan
        postBaseline = 30;  % maximum baseline (seconds) after final task period ends
        resetTime = false;  % flag indicating whether the time vector should be reset so t(1)=0 (default: false)
        Trim_Auxillary_Data=true;
    end
    
    methods
        function obj = TrimBaseline( prevJob )
           obj.name = 'Trim Pre/Post Baseline';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            
            if isempty(obj.preBaseline)
                obj.preBaseline = inf;
            end
            if isempty(obj.postBaseline)
                obj.postBaseline = inf;
            end
            
            for i = 1:numel(data)
                
                d = data(i).data;
                t = data(i).time;
                
                % get stims 
                stims = data(i).stimulus.values;
                
                % get stim vectors from stims
                o=[];
                du=[];
                for j = 1:length(stims)
                    if(isa(stims{j},'nirs.design.StimulusEvents'))
                        o=[o; stims{j}.onset(:)];
                        du=[du; stims{j}.dur(:)];
                    elseif(isa(stims{j},'nirs.design.StimulusVector'))
                        s = stims{j}.getStimVector( t );
                        onset = t(find(s~=0,1,'first'));
                        dur = t(find(s~=0,1,'last')) - onset;
                        o=[o; onset];
                        du=[du; dur];
                    else
                        error('Stimulus type not recognized: %s',class(stims{j}));
                    end
                end
                
                t_min=min(o)- obj.preBaseline;
                n=min(length(o),length(du));
                o(n+1:end)=[]; du(n+1:end)=[];
                t_max=max(o+du)+ obj.postBaseline;

                % extract data from the time inverval t_min to t_max
                lst = t >= t_min & t <= t_max;
                t = t(lst);
                d = d(lst,:);
                
                if(obj.Trim_Auxillary_Data && isa(data(i),'nirs.core.Data') && data(i).auxillary.count>0)
                    for j=1:data(i).auxillary.count
                        key= data(i).auxillary.keys{j};
                        try
                            dd=data(i).auxillary(key);
                            dd.data=dd.data(lst,:);
                            dd.time=dd.time(lst);
                            data(i).auxillary(key)=dd;
                        end
                    end
                end
                
                
                % Reset time so t(1)=0
                if obj.resetTime
                    
                    tmin = min(t);
                    t = t - tmin;
                    conds = data(i).stimulus.keys;
                    for j = 1:length(conds)
                        stims = data(i).stimulus(conds{j});
                        if(isa(stims,'nirs.design.StimulusEvents'))
                            stims.onset = stims.onset - tmin;
                            data(i).stimulus(conds{j}) = stims;
                        end
                    end
                    
                    if(obj.Trim_Auxillary_Data && isa(data(i),'nirs.core.Data') && data(i).auxillary.count>0)
                        for j=1:data(i).auxillary.count
                            key= data(i).auxillary.keys{j};
                            try
                                dd=data(i).auxillary(key);
                                dd.time=dd.time-tmin;
                                data(i).auxillary(key)=dd;
                            end
                        end
                    end
                end
                
                data(i).data = d;
                data(i).time = t;
                
            end
        end
    end
    
end

