classdef ChannelStatsROC
%% ChannelStatsROC - This class will perform ROC analysis. 
%   Takes in a function to simulate data and an analysis pipeline that ends
%   with a ChannelStats object.
    
    properties
        simfunc  = @nirs.testing.simData
        pipeline
    end
    
    properties (SetAccess = protected)
       truth
       pvals
       types
    end
    
    methods
        % constructor
        function obj = ChannelStatsROC( pipeline )
           if nargin < 1
               p = nirs.modules.Resample();
               p = nirs.modules.OpticalDensity(p);
               p = nirs.modules.BeerLambertLaw(p);
               p = nirs.modules.AR_IRLS(p);
               p.verbose = false;
               
               obj.pipeline = p;
           end
        end
        
        function obj = run(obj, iter)
            for i = 1:iter
               [data, truth] = obj.simfunc();
               
               % pipeline stats
               stats = obj.pipeline.run(data);
               
               % multivariate joint hypothesis testing
               fstats = stats.jointTest();
               
               % types
               types = unique(stats.variables.type, 'stable');
               
               t = []; p = [];
               for j = 1:length(types)
                   lst = strcmp(types(j), stats.variables.type);
                   
                   t(:,j) = truth(lst);
                   p(:,j) = stats.p(lst);
               end
               
               t(:, end+1) = sum(t, 2) > 0;
               p(:, end+1) = fstats.p;
               
               obj.truth = [obj.truth; t];
               obj.pvals = [obj.pvals; p];
               
               obj.types = [types; {'joint'}];
            end
            
            disp( ['Finished iter: ' num2str(i)] )
        end
        
        function draw(obj)
            colors = lines(length(obj.types));
            
            figure, hold on
            for i = 1:length(obj.types)
               [tp, fp, phat] = nirs.testing.roc(obj.truth(:, i), obj.pvals(:, i));
               plot(fp, tp, 'Color', colors(i,:))
            end
            xlabel('False Positive Rate')
            ylabel('True Positive Rate')
            legend(obj.types)
            
            figure, hold on
            for i = 1:length(obj.types)
               [tp, fp, phat] = nirs.testing.roc(obj.truth(:, i), obj.pvals(:, i));
               plot(phat, fp, 'Color', colors(i,:))
            end
            ylabel('False Positive Rate')
            xlabel('Estimated FPR (p-value)')
            legend(obj.types)
        end
        
        function sensitivity( obj, pval )
            
        end
        
        function specificity( obj, pval )
            
        end
        
        function auc( obj, fpr_max )
            
        end
        
    end
    
end

