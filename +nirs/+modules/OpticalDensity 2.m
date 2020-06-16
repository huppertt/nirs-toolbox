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
                

                m = mean( d, 1 );
                
                d=-log(d./(ones(size(d,1),1)*m));
                
                %d = bsxfun( @plus, -log(d), log(m) );
                
                if(~isreal(d))
                   % disp('Warning: negative intensities encountered');
                   % d= abs(d);
                   % d=max(d,eps(1));
                end
                
                data(i).data = real(d);
            end
        end
    end
    
end

