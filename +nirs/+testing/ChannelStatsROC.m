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
        dataset
        artfunc
        pipeline
    end
    properties(Hidden = true)
        beta
    end
    
    properties (SetAccess = protected)
       truth
       pvals
       types
    end
    

    methods
        % constructor
        function obj = ChannelStatsROC( pipeline, simfunc )
           if nargin < 1 || isempty(pipeline)
               p = nirs.modules.Resample();
               p = nirs.modules.OpticalDensity(p);
               p = nirs.modules.BeerLambertLaw(p);
               p = nirs.modules.GLM(p);
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
            if ~isempty(obj.dataset)
                if (~ismember('data', fieldnames(obj.dataset)))
                    error('No data feild in dataset!')
                end
                if (~ismember('truth', fieldnames(obj.dataset)))
                    error('No truth feild in dataset!')
                end
                if (iter > length(obj.dataset.data))
                    warning('Requested number of iterations is greater than data size. Discard excessive iterations.');
                    iter = length(obj.dataset.data);
                end
            end
            for idx = 1:iter
                if ~isempty(obj.dataset)
                    data = obj.dataset.data(idx);
                    truth = obj.dataset.truth(:, idx);
                else
                    [data, truth] = obj.simfunc();
                end
               
               if ~isempty(obj.artfunc)
                   data = obj.artfunc(data);
               end
               if(length(truth)>height(data(1).probe.link))
                    % Image based methods
                   error('use nirs.testing.ImageStatsROC for images');
               else
                    [~,i]=sortrows(data(1).probe.link,{'type','source','detector'});
                    truth=truth(i);
                    data=data.sorted({'type','source','detector'});
                    isimage=false;
               end
               % pipeline stats
               
               truth_save = truth;
               
               T=[];
               P=[];
               B=[];
               Types={};
               
               for i=1:length(obj.pipeline)
                   truth = truth_save;
                   if(iscell(obj.pipeline))
                       stats = obj.pipeline{i}.run(data);
                   else
                    stats = obj.pipeline(i).run(data);
                   end
                   
                   
                   if(ismember('ShortSeperation',data(1).probe.link.Properties.VariableNames))
                      if(any(stats(1).probe.link.ShortSeperation))
                           truth = truth(~data(1).probe.link.ShortSeperation);
                           job=nirs.modules.RemoveShortSeperations;
                           stats=job.run(stats);
                       end
                   end
                   
                   if ~isempty(strfind(stats(1).variables.cond,'Aux'))
                       job=nirs.modules.RemoveAuxFromStats;
                       stats=job.run(stats);
                   end
                   
                   if(length(truth)>height(data(1).probe.link))
                       stats=sorted(stats,{'cond','type','VoxID'});
                   else
                       stats=sorted(stats,{'type','source','detector'});
                   end
                   
                   % types
                   types = unique(stats.variables.type, 'stable');
                  
                   
                   
                   % multivariate joint hypothesis testing
                   fstats = stats.jointTest();
                   
                   t = []; p = []; betas=[];
                   for j = 1:length(types)
                       if(iscellstr(types(j)))
                       lst = strcmp(types(j), stats.variables.type);
                       else
                           lst = find(types(j)==stats.variables.type);
                       end
                       
                       t(:,j) = truth(lst);
                       p(:,j) = stats.p(lst);
          
                       betas(:,j)=stats.tstat(lst); %stats.beta(lst);

                   end
                   if(~iscellstr(types(1)))
                        for i=1:length(types)
                            tt{i,1}=[num2str(types(i)) 'nm'];
                        end
                        types=tt;
                   end
                   t(:, end+1) = sum(t, 2) > 0;
                   p(:, end+1) = fstats.p;
                   
                   types=[types; {'joint'}];
                   
                   if(length(obj.pipeline)>1)
                       types=arrayfun(@(x){[x{1} '-' num2str(i)]},types);
                   end
                   
                   if(isimage)
                       
                      lstp=find(t(:,1)==1);
                      lstn=find(t(:,1)==0);
                      n=min(length(lstp),length(lstn));
                      
                      kk=randperm(length(lstp));
                      lstp=lstp(kk(1:n));
                       kk=randperm(length(lstp));
                      lstn=lstn(kk(1:n));
                      t=t([lstp lstn],:);
                      p=p([lstp lstn],:);
                       
                       
                   end
                   
                   T=[T t];
                   P=[P p];
                   B=[B betas];
                   Types={Types{:} types{:}};
                   
               end
               
               obj.truth = [obj.truth; T];
               obj.pvals = [obj.pvals; P];
               obj.beta=[obj.beta; B];
               obj.types = Types;
            
               disp( ['Finished iter: ' num2str(idx)] )
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
            
           if(strcmp(type,'seperate'))
               type=unique(utype);
                seperate=true;
           else
               seperate=false;
           end
           
            colors = lines(length(obj.types));
            
           figure, hold on;
            for i = 1:length(obj.types)
                
                if(ismember(utype{i},type))
                    if(seperate)
                        hold on;
                        subplot(1,length(type),find(ismember(type,utype{i})));
                    end
                    [tp, fp, phat] = nirs.testing.roc(obj.truth(:, i), obj.pvals(:, i));
                     plot(fp, tp, 'Color', colors(i,:));
                
                end
            end
            if(~seperate)
               
                xlabel('False Positive Rate');
                ylabel('True Positive Rate');
                legend({obj.types{ismember(utype,type)}}, 'Location', 'SouthEast');
            else
                for i=1:length(type)
                    subplot(1,length(type),i);
                    title(type{i});
                    xlabel('False Positive Rate');
                    ylabel('True Positive Rate');
                    legend({obj.types{ismember(utype,type{i})}}, 'Location', 'SouthEast');
                end
                
            end
           
            figure, hold on;
            for i = 1:length(obj.types)
                if(ismember(utype{i},type))
                    if(seperate)
                        hold on;;
                        subplot(1,length(type),find(ismember(type,utype{i}))); 
                    end
                    [tp, fp, phat] = nirs.testing.roc(obj.truth(:, i), obj.pvals(:, i));
                   plot(phat, fp, 'Color', colors(i,:));
                  
                    
                end
            end
            
             if(~seperate)
                   plot([0 1],[0 1],'Color',[.8 .8 .8],'linestyle','--')
                legend({obj.types{ismember(utype,type)} 'truth'}, 'Location', 'SouthEast');
                ylabel('False Positive Rate')
                    xlabel('Estimated FPR (p-value)')
             else
                 for i=1:length(type)
                    subplot(1,length(type),i);
                    hold on;
                    title(type{i});   
                    plot([0 1],[0 1],'Color',[.8 .8 .8],'linestyle','--');
                    legend({obj.types{ismember(utype,type{i})} 'truth'}, 'Location', 'SouthEast');
                
                    ylabel('False Positive Rate')
                    xlabel('Estimated FPR (p-value)')
                 end
                 
             end
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
        
        function obj =clear(obj)
            obj.truth=[];
            obj.pvals=[];
            obj.types={};
            obj.beta=[];
        end
        
        function out = sensitivity( obj, pval )
            for i = 1:length(obj.types)
                t = obj.truth(:,i);
                p = obj.pvals(:,i);
                lst=~isnan(t);
                t=t(lst);
                p=p(lst);
                out(i,1) = sum(t(p<pval)) / sum(t);
            end
        end
        
        function out = specificity( obj, pval )
            for i = 1:length(obj.types)
                t = obj.truth(:,i);
                p = obj.pvals(:,i);
                lst=~isnan(t);
                t=t(lst);
                p=p(lst);
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
        
        function obj = comb( obj, obj1)
            obj.truth = [obj.truth; obj1.truth];
            obj.pvals = [obj.pvals; obj1.pvals];
            
        end        
    end
    
end

