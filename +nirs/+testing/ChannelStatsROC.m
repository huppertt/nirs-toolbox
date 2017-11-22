classdef ChannelStatsROC
%% ChannelStatsROC - This class will perform ROC analysis. 
% Takes in a function to simulate data and an analysis pipeline that ends
% with a ChannelStats object.
% 
% Example 1:
%     % load some data
%     raw = nirs.io.loadDirectory( '~/resting_state', {} );
%     sd  = unique([raw(1).probe.link.source raw(1).probe.link.detector], 'rows');
%     n   = size(sd,1);
% 
%     % simulation function
%     rndData = @() raw(randi(length(raw)));
%     rndChan = @() sd(rand(n,1)<0.5,:);
%     simfunc = @() nirs.testing.simData( rndData(), [], 7, rndChan() );
% 
%     % setup ROC
%     test = nirs.testing.ChannelStatsROC();
%     test.simfunc = simfunc;
% 
%     % run ROC
%     test = test.run( 10 );
% 
%     % draw ROC curves
%     test.draw()
% 
% Example 2:
%     % pipeline
%     p = nirs.modules.Resample();
%     p = nirs.modules.OpticalDensity(p);
%     p = nirs.modules.BeerLambertLaw(p);
%     p = nirs.modules.AR_IRLS(p);
%     p.verbose = false;
%     p = nirs.modules.MixedEffects(p);
%     p.formula = 'beta ~ -1 + cond';
% 
%     test = nirs.testing.ChannelStatsROC(p, @nirs.testing.simDataSet);
% 
%     test = test.run(5);
% 
% Example 3:
%     % pipeline
%     p = nirs.modules.Resample();
%     p = nirs.modules.OpticalDensity(p);
%     p = nirs.modules.BeerLambertLaw(p);
%     p = nirs.modules.AR_IRLS(p);
%     p.verbose = false;
% 
%     p = nirs.modules.MixedEffects(p);
%     p.formula = 'beta ~ -1 + cond';
%     p.dummyCoding = 'full';
% 
%     % going to use ttest to subtract cond B from A
%     pipeline.run = @(data) p.run(data).ttest([1 -1]);
% 
%     % sim function
%     randStim = @(t) nirs.testing.randStimDesign(t, 2, 7, 2);
%     simfunc = @()nirs.testing.simDataSet([], [], randStim, [5 2]', []);
% 
%     % ROC test
%     test = nirs.testing.ChannelStatsROC(pipeline, simfunc);
% 
%     test = test.run(3);

    properties
        simfunc  = @nirs.testing.simData
        artfunc
        pipeline
    end
    
    properties (SetAccess = protected)
       truth
       pvals
       types
    end
    
    methods
        % constructor
        function obj = ChannelStatsROC( pipeline, simfunc )
           if nargin < 1
               p = nirs.modules.Resample();
               p = nirs.modules.OpticalDensity(p);
               p = nirs.modules.BeerLambertLaw(p);
               p = nirs.modules.AR_IRLS(p);
               p.verbose = false;
               
               obj.pipeline = p;
           else
               obj.pipeline = pipeline;
           end
           
           if nargin > 1
               obj.simfunc = simfunc;
           end
        end
        
        function obj = run(obj, iter)
            for i = 1:iter
               [data, truth] = obj.simfunc();
               if ~isempty(obj.artfunc)
                   data = obj.artfunc(data);
               end
               if(length(truth)>height(data(1).probe.link))
                    % Image based methods
               else
               [~,i]=sortrows(data(1).probe.link,{'type','source','detector'});
               truth=truth(i);
               end
               % pipeline stats
               
               T=[];
               P=[];
               Types={};
               
               for i=1:length(obj.pipeline)
                   
                   stats = obj.pipeline(i).run(data);
                   
                    if(length(truth)>height(data(1).probe.link))
                        stats=sorted(stats,{'cond','type','VoxID'});
                    else
                        stats=sorted(stats,{'type','source','detector'});
                    end
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
                   
                    types=[types; {'joint'}];
                   
                   if(length(obj.pipeline)>1)
                       types=arrayfun(@(x){[x{1} '-' num2str(i)]},types);
                   end
                   
                   T=[T t];
                   P=[P p];
                   Types={Types{:} types{:}};
                   
               end
               
               obj.truth = [obj.truth; T];
               obj.pvals = [obj.pvals; P];
               
               obj.types = Types;
            
               it=round(length(obj.truth)/length(T));
                disp( ['Finished iter: ' num2str(it)] )
            end
        end
        
        function draw(obj,type)
            
            utype={};
            for i=1:length(obj.types);
                if(~isempty(strfind(obj.types{i},'-')))
                    utype{i}=obj.types{i}(1:max(strfind(obj.types{i},'-'))-1);
                else
                    utype{i}=obj.types{i};
                end
            end
           if(nargin<2)
                type=unique(utype);
            end
            
            colors = lines(length(obj.types));
            
            figure, hold on
            for i = 1:length(obj.types)
                if(ismember(utype{i},type))
                    [tp, fp, phat] = nirs.testing.roc(obj.truth(:, i), obj.pvals(:, i));
                    plot(fp, tp, 'Color', colors(i,:))
                end
            end
            xlabel('False Positive Rate')
            ylabel('True Positive Rate')
            legend({obj.types{ismember(utype,type)}}, 'Location', 'SouthEast')
            
            figure, hold on
            for i = 1:length(obj.types)
                if(ismember(utype{i},type))
                    [tp, fp, phat] = nirs.testing.roc(obj.truth(:, i), obj.pvals(:, i));
                    plot(phat, fp, 'Color', colors(i,:))
                end
            end
            plot([0 1],[0 1],'Color',[.8 .8 .8],'linestyle','--')
            ylabel('False Positive Rate')
            xlabel('Estimated FPR (p-value)')
            legend({obj.types{ismember(utype,type)}}, 'Location', 'SouthEast')
        end
        
        function [tp, fp, phat] = roc( obj,flag)
            if(nargin<2)
                flag=false;
            end
            pvals=obj.pvals;
            
            if(flag)
                for i=1:size(pvals,1)
                    pvals(i,1:3:end)=nirs.math.fdr(pvals(i,1:3:end));
                    pvals(i,2:3:end)=nirs.math.fdr(pvals(i,2:3:end));
                    pvals(i,3:3:end)=nirs.math.fdr(pvals(i,3:3:end));
                    
                end
            end
            
            for i = 1:length(obj.types)
               [tp{i}, fp{i}, phat{i}] = nirs.testing.roc(obj.truth(:, i),pvals(:, i));
            end
        end
        
        function out = sensitivity( obj, pval )
            for i = 1:length(obj.types)
                t = obj.truth(:,i);
                p = obj.pvals(:,i);
                out(i,1) = sum(t(p<pval)) / sum(t);
            end
        end
        
        function out = specificity( obj, pval )
            for i = 1:length(obj.types)
                t = obj.truth(:,i);
                p = obj.pvals(:,i);
                out(i,1) = 1 - sum(~t(p<pval)) / sum(~t);
            end
        end
        
        function out = auc( obj, fpr_max )
            if nargin < 2
                fpr_max = 1;
            end
            
            for i = 1:length(obj.types)
               [tp, fp, phat] = nirs.testing.roc(obj.truth(:, i), obj.pvals(:, i));
               
               lst = fp <= fpr_max;
               out(i, 1) = sum(tp(lst) .* diff([0; fp(lst)]));
            end
            
        end
        
    end
    
end

