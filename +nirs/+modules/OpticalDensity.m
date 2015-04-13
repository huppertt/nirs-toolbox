classdef OpticalDensity < nirs.modules.AbstractModule
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
                
                d = real( -log(d) );
                m = mean( d, 1 );
                d = bsxfun( @minus, d, m );
                
                data(i).data = d;
                
            end
        end
    end
    
end

