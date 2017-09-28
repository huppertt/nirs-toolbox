classdef BaselineCorrection < nirs.modules.AbstractModule
%% BaselineCorrection - Attemps a very conservative motion correction to remove DC-shifts.
% 
% Options: 
%     tune - number of standard deviations to define an outlier

    properties
        tune = 5; % number of standard deviations to define an outlier
        PCA = false;
        verbose;
        cache_dir;  % (optional) directory to cache results (unset disables caching)
        cache_rebuild;  % (optional) force rebuild of cached results (don't load previous results from cache, only save new results)
    end
    
    methods
        function obj = BaselineCorrection( prevJob )
           obj.name = 'Correct Baseline Shifts';
           obj.verbose = false;
           obj.cache_dir='';
           obj.cache_rebuild=false;
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            
            for i = 1:numel(data)
                
                % Compute data hash and load cached result if match is found
                clear hash cache_file
                if exist(obj.cache_dir,'dir')
                    hashopt.Method = 'SHA-256';
                    hash = DataHash( { data(i).data , data(i).time , obj.tune , obj.PCA } , hashopt );
                    cache_file = fullfile( obj.cache_dir , [hash '.mat'] );
                    if ~obj.cache_rebuild && exist(cache_file,'file')
                        tmp = load(cache_file,'data');
                        data(i).data = tmp.data;
                        if obj.verbose
                            disp(['Finished ' num2str(i) ' of ' num2str(length(data)) ' (cached)']);
                        end
                        continue;
                    end
                end
                
                medians = median(data(i).data);
                data(i).data = bsxfun( @minus , data(i).data , medians );
                
                if obj.PCA
                    [data(i).data,projmat] = nirs.math.pca( data(i).data );
                end
                
                for j = 1:size(data(i).data, 2)
                    y   = data(i).data(:,j);
                    Fs  = data(i).Fs;
                    t   = data(i).time;
                    
%                     % fitting an ARI model with one diff
%                     yd = diff([0; y-m]);
%                     
%                     [a, r] = nirs.math.robust_ar_fit(yd, round(4*Fs));
                    
                    [~,~,~,ymoco] = nirs.math.robust_ari1_fit(y, round(4*Fs), obj.tune);
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
                    data(i).data(:,j) = ymoco;
                    
                end
                
                if obj.PCA
                    data(i).data = data(i).data * projmat';
                end
                
                data(i).data = bsxfun(@plus,bsxfun(@minus,data(i).data,median(data(i).data)),medians);
                
                if obj.verbose
                    disp(['Finished ' num2str(i) ' of ' num2str(length(data))]);
                end
                
                if exist('cache_file','var')
                    tmp = [];
                    tmp.data = data(i).data;
                    save(cache_file,'-struct','tmp');
                end
            end
            
        end
    end
end