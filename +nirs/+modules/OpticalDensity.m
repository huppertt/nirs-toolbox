classdef OpticalDensity < nirs.functional.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    methods

        function obj = OpticalDensity( prevJob )
           obj.name = 'Optical Density';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = execute( obj, data )
            for i = 1:length(data)
                d = data(i).data;
                
                d = real( -log(d) );
                m = mean( d, 1 );
                d = bsxfun( @minus, d, m );
                
                data(i).data = d;
                
            end
        end
        
        function options = getOptions( obj )
            options = [];
        end
           
        function obj = putOptions( obj, options )
        end
        
    end
    
end

