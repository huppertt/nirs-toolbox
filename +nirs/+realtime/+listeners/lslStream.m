classdef lslStream < handle
    properties
        Fs;
        LSLdata_StreamName=[];
        LSLmarker_StreamName=[];
        data_output;
        timeout=10;
    end
    properties(Hidden=true)
        timer=[];
        LSLdata_Stream=[];
        LSLmarker_Stream=[];
        liblsl=[];
    end
    
    methods
        
        function obj=lslStream(fs)
            if(nargin==1)
                obj.Fs=fs;
            else
                obj.Fs=4;
            end
            obj.liblsl = lsl_loadlib();

        end
        
        function set.Fs(obj,Fs)
            obj.Fs=Fs;
            obj.timer=timer('Period',1/Fs,...
                    'TasksToExecute',inf,...
                    'TimerFcn',@timer_callback,'ExecutionMode','fixedRate',...
                    'UserData',obj);

            
        end
        
        function start(obj)
            if(~isempty(obj.LSLdata_StreamName))
                result = lsl_resolve_byprop(obj.liblsl,'type',obj.LSLdata_StreamName,1,obj.timeout);
                if(~isempty(result))
                    obj.LSLdata_Stream = lsl_inlet(result{1});
                    disp('Now receiving chunked data...');
                    %FIRST CHUNK PULLED IS SEEMINGLY EMPTY, THEREFORE PULL
                    %A CHUNK DURING THE INIT:
                    obj.LSLdata_Stream.pull_chunk();
                else
                    warning(['Unable to find data stream: ' obj.LSLdata_StreamName]);
                end
            end
            if(~isempty(obj.LSLmarker_StreamName))
                result = lsl_resolve_byprop(obj.liblsl,'type',obj.LSLmarker_StreamName,1,obj.timeout);
                if(~isempty(result))
                    obj.LSLmarker_Stream = lsl_inlet(result{1});
                else
                    warning(['Unable to find marker stream: ' obj.LSLmarker_StreamName]);
                end
            end
            if(isempty(obj.LSLdata_Stream) && isempty(obj.LSLmarker_Stream))
                warning('LSL streams not defined');
                return;
            end

            start(obj.timer);
        end
        function stop(obj)
            stop(obj.timer);
        end
    end
    
    
end

function timer_callback(varargin)

    pause(0.05);
    obj=varargin{1}.UserData;
    
    if(~isempty(obj.LSLdata_Stream))
        [d,t] = obj.LSLdata_Stream.pull_chunk();
        d=d(2:end,:)'; %first channel is the frame number
        t=t';
        if(~isempty(t) && ~isempty(obj.data_output))
            obj.data_output.adddata(d,t);
        end
    end
    
    if(~isempty(obj.LSLmarker_Stream))
        [d,t] = obj.LSLmarker_Stream.pull_chunk();
        if(~isempty(t) && ~isempty(obj.data_output))
            for i=1:length(t)
                 obj.data_output.addevent(['Channel_' num2str(d)],t(i));
            end
        end
    end

end