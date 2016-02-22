classdef KurtoisFilter < nirs.modules.AbstractModule
%% KurtoisFilter - Removes principal components reducing spatial covariance.
%
% Options:
%     ncomp - % number of components to remove

    methods

        function obj = KurtoisFilter( prevJob )
           obj.name = 'Remove Principal Components with high Kurtois';
           
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
                
                % first pass remove spikes
                k=kurtosis(u,[],1)';
                w = mad(k, 0) / 0.6745;
                k = k/w/4.685;
                
                w = (1 - k.^2) .* (k < 1 & k > -1);
                s=s.*w;
                d = u*diag(s)*v';
                
                % add mean back
                d = bsxfun(@plus, d, m);
                
                % put back
                data(i).data = d;                
            end
        end
    end
    
end

