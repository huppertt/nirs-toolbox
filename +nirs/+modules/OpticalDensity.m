classdef OpticalDensity < nirs.modules.AbstractModule
%% OpticalDensity - Converts raw data to optical density.
% 
% dOD = -log( raw/mean(raw) )

    methods
        function obj = OpticalDensity( prevJob )
           obj.name = 'Optical Density';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:length(data)
                d = data(i).data;
                
                m = mean( d, 1 );
                d = bsxfun( @plus, -log(d), log(m) );
                
                data(i).data = real(d);
            end
        end
    end
    
end

