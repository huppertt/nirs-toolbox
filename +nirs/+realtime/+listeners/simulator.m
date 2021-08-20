classdef simulator < handle
    properties
        Fs=4;
        datasource;
        data_output;
    end
    properties(Hidden=true)
        timer;
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
            obj.Fs=obj.datasource.Fs;
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
    obj=varargin{1}.UserData;
    idx=varargin{1}.TasksExecuted;
    if(~isempty(obj.data_output))
        obj.data_output.adddata(obj.datasource.data(idx,:),obj.datasource.time(idx));
    end
end