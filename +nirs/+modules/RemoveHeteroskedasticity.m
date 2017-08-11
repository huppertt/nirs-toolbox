classdef RemoveHeteroskedasticity < nirs.modules.AbstractModule
    %% RemoveHeteroskedasticity - Remove Heteroskedasticity in the innovations of a signal
    %
  
    methods
        function obj = RemoveHeteroskedasticity( prevJob )
            obj.name = 'Remove Heteroskedasticity';
            
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            for i = 1:length(data)
                 data(i).data=nirs.math.normrootstationarity_ar(data(i).data);
            end
        end
    end
    
end

