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
            for i = 1:numel(data)
                d = data(i).data;
                
                if(any(d(:)<=0))
                    warning('negative intensities encountered');
                    d= abs(d);
                    d=max(d,eps(1));
                end
                
                m = mean( d, 1 );
                
                
                
                d = bsxfun( @plus, -log(d), log(m) );
                
                data(i).data = real(d);
            end
        end
    end
    
end

