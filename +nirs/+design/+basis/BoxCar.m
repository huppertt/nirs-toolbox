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
            for i=1:Fs*obj.irf_dur; s(lst+i,:)=max(s(lst+i,:),s(lst,:)); end;
            
            out = [zeros(nLag,1); s(1:end-nLag)];
                        
        end
        
    end
    
end

