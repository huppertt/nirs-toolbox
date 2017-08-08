classdef PartitionedEnsemble < classreg.learning.partition.PartitionedModel
%PartitionedEnsemble Cross-validated model.
%   PartitionedEnsemble is the super class for cross-validated ensembles.

%   Copyright 2010-2013 The MathWorks, Inc.

    
    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        %TRAINABLE Full models trained on cross-validation folds.
        %   The Trainable property is a cell array of ensembles trained on
        %   cross-validation folds. Every ensemble is full, that is, it stores X, Y
        %   and W data it has been trained on.
        %
        %   See also ClassificationPartitionedEnsemble, RegressionPartitionedEnsemble.
        Trainable;
        
        %NUMTRAINEDPERFOLD Number of trained learners per cross-validation fold.
        %   The NumTrainedPerFold property is a vector with KFold elements. Every
        %   element shows the number of trained weak learners in this fold.
        %
        %   See also ClassificationPartitionedEnsemble,
        %   RegressionPartitionedEnsemble, PartitionedModel/KFold.
        NumTrainedPerFold;
    end
    
    properties(GetAccess=public,SetAccess=protected,Hidden=true,Dependent=true)
        NTrainedPerFold;
    end
    
    properties(GetAccess=public,SetAccess=protected,Dependent=true,Hidden=true)
        Combiner; % Cell array of combiners across all folds
    end
       
    methods
        function trainable = get.Trainable(this)
            trainable = this.Ensemble.Trainable;
        end
        
        function tfold = get.NumTrainedPerFold(this)
            tfold = this.NTrainedPerFold;
        end
        
        function tfold = get.NTrainedPerFold(this)
            kfold = length(this.Ensemble.Trained);
            tfold = zeros(1,kfold);
            for k=1:kfold
                tfold(k) = length(this.Ensemble.Trained{k}.Trained);
            end
        end
        
        function comb = get.Combiner(this)
            kfold = length(this.Ensemble.Trained);
            comb = cell(kfold,1);
            for k=1:kfold
                comb{k} = this.Ensemble.Trained{k}.Impl.Combiner;
            end
        end
    end
    
    methods(Access=protected)
        function this = PartitionedEnsemble()
            this = this@classreg.learning.partition.PartitionedModel();
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.partition.PartitionedModel(this,s);
            s.NumTrainedPerFold = this.NumTrainedPerFold;
        end

        function trainable = resumePartitionedWithPrint(this,nlearn,nprint)
            trainable = this.Ensemble.Trainable;
            if isempty(trainable)
                error(message('stats:classreg:learning:partition:PartitionedEnsemble:resumePartitionedWithPrint:NoTrainableLearners'));
            end
            for t=1:numel(trainable)
                trainable{t} = resume(trainable{t},nlearn);
                if mod(t,nprint)==0
                    fprintf(1,'%s%i\n',this.Ensemble.ModelParams.PrintMsg,t);
                end
            end
        end
        
        function [ensembleMode,folds,extraArgs] = checkEnsembleFoldArgs(this,varargin)
            kfold = length(this.Ensemble.Trained);
            
            args = {'mode'     'folds' 'learners'};
            defs = {'average' 1:kfold          []};
            [cvmode,folds,learners,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            if ~isempty(learners)
                error(message('stats:classreg:learning:partition:PartitionedEnsemble:checkEnsembleFoldArgs:LearnersNoop'));                
            end
            
            cvmode = lower(cvmode);
            if     strncmpi(cvmode,'average',length(cvmode))
                ensembleMode = 'ensemble';
            elseif strncmpi(cvmode,'individual',length(cvmode))
                ensembleMode = 'individual';
            elseif strncmpi(cvmode,'cumulative',length(cvmode))
                ensembleMode = 'cumulative';
            else
                error(message('stats:classreg:learning:partition:PartitionedEnsemble:checkEnsembleFoldArgs:BadMode'));
            end
            
            if islogical(folds)
                if ~isvector(folds) || length(folds)~=kfold
                    error(message('stats:classreg:learning:partition:PartitionedEnsemble:checkEnsembleFoldArgs:BadLogicalIndices', kfold));
                end
                folds = find(folds);
            end
            if ~isnumeric(folds) || ~isvector(folds) || min(folds)<=0 || max(folds)>kfold
                error(message('stats:classreg:learning:partition:PartitionedEnsemble:checkEnsembleFoldArgs:BadNumericIndices', kfold));
            end
            folds = ceil(folds);
            for k=1:length(folds)
                if this.Ensemble.Trained{folds(k)}.NTrained==0
                    warning(message('stats:classreg:learning:partition:PartitionedEnsemble:checkEnsembleFoldArgs:EmptyFolds'));
                end
            end
        end
    end

    methods(Static,Hidden)
        function [combiner,score] = predictKfoldWithCache(combiner,X,...
                t,useNfort,useDfort,trained,classnames,nonzeroProbClasses,...
                defaultScore)
            % Get number of partitions
            kfold = length(trained);

            % Get number of classes
            K = length(classnames);
            doclass = true;
            if K==0
                doclass = false;
                K = 1;
            end
            
            % Init scores
            N = size(X,1);
            score = NaN(N,K);
            
            % Loop over partitions
            for k=1:kfold
                useNfortK = useNfort(:,k); % observations for learner t
                useDfortTK = useDfort(:,t,k); % predictors for learner t
                weak = trained{k}.Impl.Trained{t};

                % If not all predictors have been used to train learner t,
                % treat this as proof that the ensemble has been grown by
                % random subspace. In that case, do not compute predictions
                % for observations with missing inputs. To compute
                % prediction for an observation with missing values,
                % subspace ensemble averages over learners trained on
                % non-missing inputs. This is consistent with what
                % CompactEnsemble does.
                if ~all(useDfortTK)
                    goodobs = ~any(isnan(X(useNfortK,useDfortTK)),2);
                    idxNfortK = find(useNfortK);
                    idxobs = idxNfortK(goodobs);
                    goodUseNfortK = false(N,1);
                    goodUseNfortK(idxobs) = true;
                else
                    goodobs = true(sum(useNfortK),1);
                    goodUseNfortK = useNfortK;
                end
                
                if doclass
                    [~,pos] = ismember(weak.ClassSummary.ClassNames,classnames);
                    [~,score(goodUseNfortK,pos)] = predict(weak,X(goodUseNfortK,useDfortTK));
                else
                    score(goodUseNfortK) = predict(weak,X(goodUseNfortK,useDfortTK));
                end
                
                combiner{k} = updateCache(combiner{k},score(useNfortK,:),t,goodobs);
                score(useNfortK,:) = cachedScore(combiner{k});
            end
            
            % Assign scores only for classes with non-zero probability
            if doclass
                tf = ismember(classnames,nonzeroProbClasses);
                score(:,~tf) = defaultScore;
            end
        end
                        
        function tf = usePredInFold(folds,T,D,trained)
            Kfold = numel(folds);
            tf = true(D,T,Kfold);
            for k=1:Kfold
                if ~isempty(trained{k}.UsePredForLearner)
                    tf(:,:,k) = trained{k}.UsePredForLearner(:,1:T);
                end
            end
        end
    end
    
end
