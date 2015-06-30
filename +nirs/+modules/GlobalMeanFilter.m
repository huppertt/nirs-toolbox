classdef GlobalMeanFilter < nirs.modules.AbstractModule
    
    methods

        function obj = GlobalMeanFilter( prevJob )
           obj.name = 'Remove Global Mean';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:length(data)
                
                % get data
                d = data(i).data;
                
                types = data(i).probe.link.type;

                if ~iscellstr(types)
                    types = arrayfun(@(x){num2str(x)}, types);
                end
                
                utypes  = unique(types);
                            
                for j = 1:length( utypes )
                    lst = strcmp(utypes(j), types);
                    
                    % global mean
                    s = mad(d(:,lst),0,1)' / 0.06745;
                    m = (1./s) \ (diag(1./s)*d(:,lst)');
                    m = m' - mean(m);
                    
                    % regressor
                    X = [m.^0 m.^1 m.^2 m.^3];
                    
                    % removal
                    tmp(:,lst) = d(:,lst) - X * (X\d(:,lst));
                end
                
                % put back
                data(i).data = d;
                
            end
        end
    end
    
end

