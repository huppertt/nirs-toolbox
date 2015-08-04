classdef BaselineCorrection < nirs.modules.AbstractModule
%% BaselineCorrection - Attemps a very conservative motion correction to remove DC-shifts.
% 
% Options: 
%     tune - number of standard deviations to define an outlier

    properties
        tune = 5; % number of standard deviations to define an outlier
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
                    
                    m = median(y);
                    
%                     % fitting an ARI model with one diff
%                     yd = diff([0; y-m]);
%                     
%                     [a, r] = nirs.math.robust_ar_fit(yd, round(4*Fs));
                    
                    [~,~,~,ymoco] = nirs.math.robust_ari1_fit(y-m, round(4*Fs), obj.tune);
%                     % weights
%                     s = mad(r, 0) / 0.6745;
%                     w = abs(r/s/obj.tune) < 1;
%                     
%                     % filter
%                     f = [1; -a(2:end)];                  
%                     
%                     ymoco = cumsum( filter(1, f, w.*r) );
%                     
%                     % remove low poly behaviour
%                     X = nirs.design.trend.legendre(t, 2);
%                     ymoco = ymoco - X*(X\ymoco);
%                     
                    ymoco = ymoco + m;
                    
                    data(i).data(:,j) = ymoco;
                end
            end
            
        end
    end
end