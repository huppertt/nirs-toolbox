classdef BaselineCorrection < nirs.modules.AbstractModule
    
    properties
        tune = 5;
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
                    
%                     [a, r] = nirs.math.ar_fit(diff([0; y]), round(4*Fs));
%                     
%                     lst = abs(r)./(mad(r,0)/0.6745) > obj.tune;
%                     r(lst) = 0;
%                     
%                     f = [1; -a(2:end)];
%                     ymoco2 = cumsum(filter(1, f, r));
                    
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