classdef BandPass < handle
    % real-time implementation of a band-pass filter
    
    properties
       lowpass;
       highpass;
       fs;
    end
    
    properties(Hidden=true);
       model=40;
       previousdata=[];
    end
    properties(Dependent=true)
         fb;
    end
    
    methods
        function obj=BandPass(lpf,hpf,fs)
            if(nargin>0)
                obj.lowpass=lpf;
            end
            if(nargin>1)
                obj.highpass=hpf;
            end
            if(nargin>2)
                obj.fs=fs;
            end
            
        end
        function fb =get.fb(obj)
            if(isempty(obj.fs))
                fb=[];
                return;
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
        
        function dOut = update(obj,d,t)
            dOut=zeros(size(d));
            for i=1:size(d,1)
                if(isempty(obj.previousdata))
                    obj.previousdata=zeros(obj.model,size(d,2));
                    obj.previousdata=ones(obj.model,1)*d(1,:);
                end
                obj.previousdata(2:end)=obj.previousdata(1:end-1);
                obj.previousdata(1,:)=d(i,:);
                dOut(i,:)=sum((obj.fb'*ones(1,size(d,2))).*obj.previousdata,1)';
            end
        end
            
    end
    
end