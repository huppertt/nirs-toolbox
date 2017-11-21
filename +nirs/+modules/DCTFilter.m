classdef DCTFilter < nirs.modules.AbstractModule
%%  DCT Filter - Apply DCT filter to data 
%
% 
% Options: 
%       cutoff = 128; % Period in seconds of filter cutoff
%     

    properties
        cutoff = 128;
    end
    
    methods

        function obj =  DCTFilter( prevJob )
           obj.name = 'Apply DCT filter';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            

            for sub = 1:length(data)
                
                X0 = nirs.design.trend.dctmtx( data(sub).time , 1/obj.cutoff );
                X0 = X0(:,2:end);

                data(sub).data = data(sub).data - X0 * (X0 \ data(sub).data);

            end
        end
    end
end
