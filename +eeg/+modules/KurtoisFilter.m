classdef KurtoisFilter < nirs.modules.AbstractModule
%% KurtoisFilter - Removes principal components reducing spatial covariance.
%
% Options:
%     ncomp - % number of components to remove
    properties
        model = 'PCA'; % use ICA or PCA
    end
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
                
                if(strcmp(obj.model,'PCA'))
                    % svd
                    [u, s, v] = svd(d,'econ');
                    s = diag(s);
                elseif(strcmp(obj.model,'ICA'))
                    [u,v,~]=fastica(d','approach', 'symm', 'g', 'tanh', 'verbose', 'off', 'displayMode', 'off','stabilization','on');
                    u=u';
                    s=ones(size(u,2),1);
                else
                    error('unrecognized model [ICA or PCA]');
                end
                % first pass remove spikes
                k=kurtosis(u,0,1)';
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

