classdef Vestibular
% Canonical HRF response
    
    properties
        peakTime    = 4;
        uShootTime  = 16;
        peakDisp    = 1;
        uShootDisp  = 1;
        ratio       = 1/6; 
        duration    = 32;
        elongate    = 40;
       
        incDeriv = false;
    end
    
    methods
        function out = convert( obj, s, t )
        
            % params
            a1 = obj.peakTime;
            a2 = obj.uShootTime;
            b1 = obj.peakDisp;
            b2 = obj.uShootDisp;
            c  = obj.ratio;
            elongate = obj.elongate;
            % sampling freq
            Fs = 1/(t(2)-t(1));

            % time vector
            t = (0:1/Fs:[obj.duration+elongate])';

            % impulse response
            h = obj.getImpulseResponse( a1, b1, a2, b2, c, t,elongate*Fs );
            
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
        function h = getImpulseResponse( a1, b1, a2, b2, c, t,elongate )
            h = b1^a1*t.^(a1-1).*exp(-b1*t)/gamma(a1) - c*b2^a2*t.^(a2-1).*exp(-b2*t)/gamma(a2);
            h=conv(h,ones(elongate,1));
            h=h(1:length(t));
            h = h / sum(h);
        end
    end
            
end