classdef KalmanGLM < nirs.realtime.modules.AbstractModule
   properties
       Pmax=8;
       Q=0;
       returntype='tstat'; % or beta
       number_conditions=1;
   end
   properties(Hidden=true)
        kfarirls; 
        Fs=1;
        lastT=[];
   end
    
    methods
        function obj = KalmanGLM(prevJob)
            obj.name='RT-GM';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function obj=resetThis(obj)
            obj.kfarirls=[];
        end
        
        function [d,t,probe,stimulus] = updateThis(obj,d,t,probe,stimulus)
            try
            nchan=height(probe.link);
            if(isempty(obj.kfarirls))
                obj.lastT=t;
                for i=1:nchan
                    modelkf=nirs.realtime.util.RobustKalmanFilter(obj.Q);
                    arkf=nirs.realtime.util.KalmanAR(obj.Pmax);
                    obj.kfarirls{i,1}=nirs.realtime.util.KalmanARWLS(modelkf,arkf);
                end
            else
                obj.Fs=1./(t-obj.lastT);
                obj.lastT=t;
            end


            if(isa(stimulus,'Dictionary'))
                if(isempty(stimulus))
                    X=zeros(1,number_conditions+1);
                    X(end)=1;
                else
                    X=nirs.design.createDesignMatrix(stimulus,[t-1/obj.Fs t]);
                    if(size(X,2)<obj.number_conditions+1)
                        X(:,end+1:obj.number_conditions+1)=0;
                    end
                    X(:,end)=1;
                    X=X(end,:);
                end
            else
                X=stimulus;
            end
            for i=1:nchan
                obj.kfarirls{i}.update(d(i),X);
                if(strcmp(obj.returntype,'tstat'))
                    d(i)=obj.kfarirls{i}.tstat();
                else
                    d(i)=obj.kfarirls{i}.B();
                end
            end
            end

            
        end
        
    end
    
    
end