classdef BoxCar
% Boxcar (1 or 0) basis with optional lag time.
    
    properties
        lagTime = 0;
    end
    
    methods
        function out = convert( obj, s, t )
            
            Fs = 1/(t(2)-t(1));
            nLag = Fs*obj.lagTime;
            
            out = [zeros(nLag,1); s(1:end-nLag)];
                        
        end
        
    end
    
end

