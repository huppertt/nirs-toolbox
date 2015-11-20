classdef OpticalDensity2Intensity < nirs.modules.AbstractModule
%% OpticalDensity - Converts optical density back to raw data.
% 
% dOD = -log( raw/mean(raw) )

    methods
        function obj = OpticalDensity2Intensity( prevJob )
           obj.name = 'dOD to Intensity';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:length(data)
                d = data(i).data;
                
                m=1./mad(d);
                d = exp(-d)*diag(m);
                
                data(i).data = real(d);
            end
        end
    end
    
end

