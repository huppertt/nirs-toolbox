classdef nonlinearHRF
% nonlinear Canonical HRF response
    
    properties
        peakTime    = struct('value',4,'fit',true,'lower',1,'upper',6);
        uShootTime  = struct('value',16,'fit',true,'lower',6,'upper',20);
        peakDisp    = struct('value',1,'fit',true,'lower',.5,'upper',6);
        uShootDisp  = struct('value',1,'fit',true,'lower',.5,'upper',6);
        ratio       = struct('value',1/6,'fit',true,'lower',.1,'upper',.5); 
        duration    = 32;
        
        incDeriv = false;
    end
    
    methods
        function x = invtform_params(obj)
            
            flds={'peakTime','uShootTime','peakDisp','uShootDisp','ratio'};
            cnt=1;
            for idx=1:length(flds)
                if(obj.(flds{idx}).fit)
                    x(cnt,1)=obj.unboundparam(obj.(flds{idx}).value,...
                        obj.(flds{idx}).lower,obj.(flds{idx}).upper);
                    cnt=cnt+1;
                end
            end
            
        end
        function obj = fwdtform_params(obj,x)
            flds={'peakTime','uShootTime','peakDisp','uShootDisp','ratio'};
            cnt=1;
            for idx=1:length(flds)
                if(obj.(flds{idx}).fit)
                    obj.(flds{idx}).value=obj.boundparam(x(cnt),...
                        obj.(flds{idx}).lower,obj.(flds{idx}).upper);
                    cnt=cnt+1;
                end
            end
        end
        
        
        function out = convert( obj, s, t )
        
            % params
            a1 = obj.peakTime.value;
            a2 = obj.uShootTime.value;
            b1 = obj.peakDisp.value;
            b2 = obj.uShootDisp.value;
            c  = obj.ratio.value;
            
            % sampling freq
            Fs = 1/(t(2)-t(1));

            % time vector
            t = (0:1/Fs:obj.duration)';

            % impulse response
            h = obj.getImpulseResponse( a1, b1, a2, b2, c, t );
            
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
