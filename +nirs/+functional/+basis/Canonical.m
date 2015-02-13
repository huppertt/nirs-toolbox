classdef Canonical
    %UNTITLED7 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        peakTime    = 6;
        uShootTime  = 16;
        peakDisp    = 1;
        uShootDisp  = 1;
        ratio       = 1/6; 
        duration    = 32;
    end
    
    methods
        function out = convert( obj, s, t )
        
            a1 = obj.peakTime;
            a2 = obj.uShootTime;
            b1 = obj.peakDisp;
            b2 = obj.uShootDisp;
            c = obj.ratio;
            
            Fs = 1/(t(2)-t(1));

            t = (0:1/Fs:obj.duration)';

            h = b1^a1*t.^(a1-1).*exp(-b1*t)/gamma(a1) - c*b2^a2*t.^(a2-1).*exp(-b2*t)/gamma(a2);
            h = h / sum(h);
            
            s1 = s(1,:);
            s = bsxfun( @minus, s, s1 );
            s = filter( h, 1, s );
            
            s = bsxfun( @plus, s, s1 );
            
            out = s;
            
        end
        
    end
    
end

