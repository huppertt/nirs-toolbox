classdef BandPass < nirs.realtime.modules.AbstractModule
    % real-time implementation of a band-pass filter
    
    properties
       lowpass=.5;
       highpass=0.016;
      
    end
    
    properties(Hidden=true);
       model=41;
       previousdata=[];
       fs=[];
       firstsample=[];
    end
    properties(Dependent=true, Hidden=true)
         fb=[];
        
    end
    
    methods
        function obj=BandPass(prevJob)
            obj.name='RT-BandPass Filter';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function obj=resetThis(obj)
            obj.previousdata=[];
            obj.fs=[];
            obj.firstsample=[];
        end
        
        function fb =get.fb(obj)
            if(isempty(obj.fs) || obj.fs==inf || obj.fs<=0)
                fb=[];
                return;
            end
            
            if(~isempty(obj.lowpass) && obj.lowpass>=obj.fs/2)
                obj.lowpass=[];
                warning('lowpass filter exceeds Nyquist');
            end
            if(~isempty(obj.highpass) && obj.highpass>=obj.fs/2)
                obj.highpass=[];
                warning('highpass filter exceeds Nyquist');
            end
            
            if(isempty(obj.lowpass) & ~isempty(obj.highpass))
                fb=fir1(obj.model-1,obj.highpass*2/obj.fs,'high');
            elseif(~isempty(obj.lowpass) & isempty(obj.highpass))
                fb=fir1(obj.model-1,obj.lowpass*2/obj.fs,'low');
            elseif(~isempty(obj.lowpass) & ~isempty(obj.highpass))
                fb=fir1(obj.model-1,[obj.highpass obj.lowpass]*2/obj.fs);
            else
                fb=[];
            end
        end
        
        function [dOut,t,probe,stimulus] = updateThis(obj,d,t,probe,stimulus)
            
            if(isempty(obj.firstsample) & isempty(obj.fs))
                obj.firstsample=t;
                dOut=d; 
                return;
            end
            if(isempty(obj.fs))
                obj.fs=1/(t-obj.firstsample);
            end
            
            if(isempty(obj.fb))
                dOut=d; 
                return;
            end
            
            dOut=zeros(size(d));
            for i=1:size(d,1)
                if(isempty(obj.previousdata))
                    obj.previousdata=zeros(obj.model,size(d,2));
                    obj.previousdata=ones(obj.model,1)*d(1,:);
                end
                obj.previousdata(2:end,:)=obj.previousdata(1:end-1,:);
                obj.previousdata(1,:)=d(i,:);
                dOut(i,:)=sum((obj.fb'*ones(1,size(d,2))).*obj.previousdata,1)';
            end
        end
            
    end
    
end