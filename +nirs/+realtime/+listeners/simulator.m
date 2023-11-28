classdef simulator < handle
    properties
        Fs=4;
        datasource;
        data_output;
    end
    properties(Hidden=true)
        timer;
        stimevents;
    end
    
    methods
        
        function obj=simulator(data)
            if(nargin==1)
                obj.datasource=data;
            end
            
        end
        
        function set.Fs(obj,Fs)
            obj.Fs=Fs;
            try
                obj.timer=timer('Period',1/Fs,...
                    'TasksToExecute',length(obj.datasource.time),...
                    'TimerFcn',@timer_callback,'ExecutionMode','fixedRate',...
                    'UserData',obj);
            end
            
        end
        
        function set.datasource(obj,data)
            if(isa(data,'nirs.core.Data'))
                obj.datasource=data;
            else
                obj.datasource=feval(data);
            end
            obj.Fs=floor(obj.datasource.Fs);

            onsets=cell(0,2);
            offsets=cell(0,2);
            stim=obj.datasource.stimulus;
            for i=1:stim.count;
                ss=stim(stim.keys{i});
                for j=1:length(ss.onset)
                    onsets{end+1,1}=ss.onset(j);
                    onsets{end,2}=stim.keys{i};
                    offsets{end+1,1}=ss.onset(j)+ss.dur(j);
                    offsets{end,2}=stim.keys{i};
                end
            end
             obj.stimevents.onsets=onsets;
             obj.stimevents.offsets=offsets;
             
        end
        
        function start(obj)
           obj.data_output.probe=obj.datasource.probe;
           start(obj.timer);
        end
        function stop(obj)
            stop(obj.timer);
            
        end
    end
    
    
end

function timer_callback(varargin)
try
    obj=varargin{1}.UserData;
    idx=varargin{1}.TasksExecuted;
    if(~isempty(obj.data_output))
        obj.data_output.adddata(obj.datasource.data(idx,:),obj.datasource.time(idx));
        lstOn=find((vertcat(obj.stimevents.onsets{:,1})>=obj.datasource.time(idx)-1/obj.Fs) & (vertcat(obj.stimevents.onsets{:,1})<obj.datasource.time(idx)+1/obj.Fs));
        lstOff=find((vertcat(obj.stimevents.offsets{:,1})>=obj.datasource.time(idx)-1/obj.Fs) & (vertcat(obj.stimevents.offsets{:,1})<obj.datasource.time(idx)+1/obj.Fs));
        
        if(~isempty(lstOn) & ~isempty(lstOff))
            obj.data_output.addevent(obj.stimevents.onsets{lstOn,2},obj.datasource.time(idx),1/obj.Fs,1);
        elseif(~isempty(lstOn))
            disp(lstOn);
            obj.data_output.addeventStart(obj.stimevents.onsets{lstOn,2},obj.datasource.time(idx),1);
        elseif(~isempty(lstOff))
            obj.data_output.addeventEnd(obj.stimevents.offsets{lstOff,2},obj.datasource.time(idx),1);
        end

    end
end
end