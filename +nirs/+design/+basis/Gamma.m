classdef Gamma
% Canonical HRF response
    
    properties
        peakTime    = 6;
        peakDisp    = 1;
        duration    = 32;
    end
    
    methods
        function obj = Gamma( ttp, dsp, dur )
            if nargin > 0, obj.peakTime = ttp; end
            if nargin > 1, obj.peakDisp = dsp; end
            if nargin > 2, obj.duration = dur; end
        end
        
        function out = convert( obj, s, t )
        
            % params
            a = obj.peakTime;
            b = obj.peakDisp;
            
            % sampling freq
            Fs = 1/(t(2)-t(1));

            % time vector
            t = (0:1/Fs:obj.duration)';

            % impulse response
            h = b^a*t.^(a-1).*exp(-b*t)/gamma(a);
            h = h / sum(h);
            
            % convert stim vectors
%             s1 = s(1,:);
%             s = bsxfun( @minus, s, s1 );
            s = filter( h, 1, s );
            
           % s = bsxfun( @plus, s, s1 );
            
            out = s;
            
        end
        
    end
    
end

