classdef BoxCar
% Boxcar (1 or 0) basis with optional lag time.
    
    properties
        lagTime = 3;
        irf_dur=5;
    end
    
    methods
        function out = convert( obj, s, t )
            
            Fs = 1/(t(2)-t(1));
            nLag = Fs*obj.lagTime;
            
            lst=find(s~=0);
            for i=1:Fs*obj.irf_dur; 
                if(max(lst+i)<length(s))
                    s(lst+i,:)=max(s(lst+i,:),s(lst,:)); 
                end
            end;
            
            out = [zeros(nLag,1); s(1:end-nLag)];
                        
        end
        function h = draw(obj)
           Fs =4; 
           t = (0:1/Fs:(obj.irf_dur+obj.lagTime+5))';
           s=zeros(size(t));
           s(1)=1;
           out = obj.convert(s, t );
           h=plot(t,out,'k');
            
            
        end
    end
    
end

