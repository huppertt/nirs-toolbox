classdef ERP
% Canonical ERP response
    
    properties
        peakTime    = .4;
        uShootTime  = 1.6;
        peakDisp    = .1;
        uShootDisp  = .1;
        ratio       = 1/4; 
        duration    = 3.2;
        
        incDeriv = false;
    end
    
    methods
        function out = convert( obj, s, t )
        
            % params
            a1 = obj.peakTime*10;
            a2 = obj.uShootTime*10;
            b1 = obj.peakDisp*10;
            b2 = obj.uShootDisp*10;
            c  = obj.ratio;
            
            % sampling freq
            Fs = 1/(t(2)-t(1));

            % time vector
            t = (0:1/Fs:obj.duration)';

            % impulse response
            h = obj.getImpulseResponse( a1, b1, a2, b2, c, t*10 );
            
%             % stupid filtering function
%             f = @(s) bsxfun( @minus, s, s(1,:) );
%             f = @(h, s) filter(h, 1, f(s));
%             f = @(h, s) bsxfun( @plus, f(h,s), s(1,:) );
            
            % convert stim vectors
            out = filter(h, 1, s);
            
            % derivatives
            if obj.incDeriv
                da = 1e-6 * a1;
                db = 1e-6 * b1;
                
                ha = (h - obj.getImpulseResponse(a1+da, b1, a2, b2, c, t) )/da;
                hb = (h - obj.getImpulseResponse(a1, b1+db, a2, b2, c, t) )/da;
                
                out = [out filter(ha, 1, s) filter(hb, 1, s)];
            end
            
        end
    end
    
    methods ( Static )
        function h = getImpulseResponse( a1, b1, a2, b2, c, t )
            h = b1^a1*t.^(a1-1).*exp(-b1*t)/gamma(a1) - c*b2^a2*t.^(a2-1).*exp(-b2*t)/gamma(a2);
            h = h / sum(h);
        end
    end
            
end

