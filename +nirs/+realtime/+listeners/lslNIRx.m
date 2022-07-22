classdef lslNIRx < handle
    properties
       data_output;
       data_type={'HbO','Hb'};
    end
    properties(Hidden=true)
        Fs;
        timer=[];
        LSLdata_Stream=[];
        LSLmarker_Stream=[];
        liblsl=[];  
        LSLdata_StreamName='Aurora';
        LSLmarker_StreamName=[]; 
        timeout=10;
        datalist;
    end
    
    methods
        
        function obj=lslNIRx(fs)
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
                result = lsl_resolve_byprop(obj.liblsl,'name',obj.LSLdata_StreamName,1,obj.timeout);
                if(~isempty(result))
                    obj.LSLdata_Stream = lsl_inlet(result{1});
                    disp('Now receiving chunked data...');
                    %FIRST CHUNK PULLED IS SEEMINGLY EMPTY, THEREFORE PULL
                    %A CHUNK DURING THE INIT:
                    obj.LSLdata_Stream.pull_chunk();
                    [d,t]=obj.LSLdata_Stream.pull_chunk();
                    n=(size(d,1)-1)/4;
                    obj.datalist=[];
                    type={};
                    link=obj.data_output.probe.link;
                    c=0;
                    if(ismember({'760'},obj.data_type))
                          obj.datalist=[obj.datalist 1:n];
                          type=[type; repmat({'760'},n,1)];
                          c=c+1;
                    end
                    if(ismember({'850'},obj.data_type))
                          obj.datalist=[obj.datalist n+1:2*n];
                            type=[type; repmat({'850'},n,1)];
                            c=c+1;
                    end
                     if(ismember({'HbO'},obj.data_type))
                          obj.datalist=[obj.datalist 2*n+1:3*n];
                            type=[type; repmat({'HbO'},n,1)];
                            c=c+1;
                     end
                     if(ismember({'Hb'},obj.data_type))
                          obj.datalist=[obj.datalist 3*n+1:4*n];
                          type=[type; repmat({'Hb'},n,1)];
                          c=c+1;
                     end
                     link=repmat(link,c,1);
                     link.type=type;
                     obj.data_output.probe.link=link;
                     
                     
                else
                    warning(['Unable to find data stream: ' obj.LSLdata_StreamName]);
                end
            end
            if(~isempty(obj.LSLmarker_StreamName))
                result = lsl_resolve_byprop(obj.liblsl,'name',obj.LSLmarker_StreamName,1,obj.timeout);
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
        d=d(1+obj.datalist,:)'; %first channel is the frame number
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