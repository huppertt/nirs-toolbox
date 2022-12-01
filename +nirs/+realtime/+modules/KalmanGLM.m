classdef KalmanGLM < nirs.realtime.modules.AbstractModule
   properties
       Pmax=8;
       Q=0;
       returntype='tstat'; % or beta
   end
   properties(Hidden=true)
        kfarirls;  
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
            
            nchan=height(probe.link);
            if(isempty(obj.kfarirls))
                for i=1:nchan
                    modelkf=nirs.realtime.util.RobustKalmanFilter(obj.Q);
                    arkf=nirs.realtime.util.KalmanAR(obj.Pmax);
                    obj.kfarirls{i,1}=nirs.realtime.util.KalmanARWLS(modelkf,arkf);
                end
                
            end
            if(isa(stimulus,'Dictionary'))
                if(isempty(stimulus))
                    X=1;
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