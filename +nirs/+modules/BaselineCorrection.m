classdef BaselineCorrection < nirs.modules.AbstractModule
    
    properties
        tune = 4.685;
    end
    
    methods
        function obj = BaselineCorrection( prevJob )
           obj.name = 'Correct Baseline Shifts';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            
            for i = 1:length(data)
                for j = 1:size(data(i).data, 2)
                    y   = data(i).data(:,j);
                    Fs  = data(i).Fs;
                    t   = data(i).time;
                    
                    [~, ~, ~, ymoco] = ...
                        nirs.math.robust_ari1_fit(y, round(4*Fs), obj.tune);
                    
                    X = nirs.design.trend.legendre(t, 3);
                    ymoco = ymoco - X*(X\ymoco);
                    
                    data(i).data(:,j) = ymoco;
                end
            end
            
        end
    end
end