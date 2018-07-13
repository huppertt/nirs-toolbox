classdef VarianceNormalize < nirs.modules.AbstractModule
%% Normalizes the variance per channel
%

    methods

        function obj = VarianceNormalize( prevJob )
           obj.name = 'Variance Normalize';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:length(data)
                
                % resample data
                d = data(i).data;
                d = d./(ones(size(d,1),1)*sqrt(var(d,[],1)));
                
                % put back
                data(i).data = d;                
            end
        end
    end
    
end

