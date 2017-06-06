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
                
                k  = size(data(sub).data, 1); 
                RT = 1/data(sub).Fs; 
                n  = fix(2*(k*RT)/obj.cutoff + 1);
                X0 = spm_dctmtx(k,n);
                X0 = X0(:,2:end);

                data(sub).data = data(sub).data - X0 * (X0' * data(sub).data);

            end
        end
    end
end
