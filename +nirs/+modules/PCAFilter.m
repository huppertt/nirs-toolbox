classdef PCAFilter < nirs.modules.AbstractModule
%% PCAFilter - Removes principal components reducing spatial covariance.
%
% Options:
%     ncomp - % number of components to remove

    properties
        ncomp = 1; % number of components to remove
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
                
                if(obj.ncomp>0 & obj.ncomp<1)
                    % if fraction, then remove that % of the noise
                    n=max(cumsum(s)/sum(s)<=obj.ncomp);
                    disp(['Removing first ' num2str(n) ' components (' num2str(cumsum(s(1:n))/sum(s)*100) '%)']);
                    
                    if(n==0)
                        disp(['     lowest component = ' num2str(cumsum(s(1))/sum(s)*100) '%)']);
                    end
                else
                    n=obj.ncomp;
                end
                
                % remove n components
                s(1:n) = 0;
                d = u*diag(s)*v';
                
                % add mean back
                d = bsxfun(@plus, d, m);
                
                % put back
                data(i).data = d;                
            end
        end
    end
    
end

