classdef PCAFilter < nirs.modules.AbstractModule
  
    properties
        ncomp = 1;
    end
    
    methods

        function obj = PCAFilter( prevJob )
           obj.name = 'Remove Principal Components';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:length(data)
                
                % resample data
                d = data(i).data;
                
                % remove mean
                m = mean(d,1);
                d = bsxfun(@minus, d, m);
                
                % svd
                [u, s, v] = svd(d,'econ');
                s = diag(s);
                
                % remove n components
                s(1:obj.ncomp) = 0;
                d = u*diag(s)*v';
                
                % add mean back
                d = bsxfun(@plus, d, m);
                
                % put back
                data(i).data = d;                
            end
        end
    end
    
end

