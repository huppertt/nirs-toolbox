classdef CompactEnsemble
%CompactEnsemble Compact ensemble.
%   CompactEnsemble is the super class for compact ensemble models.
    
%   Copyright 2010-2013 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=public,Hidden=true,Abstract=true)
        Impl;
    end

    properties(GetAccess=public,SetAccess=protected)
        %USEPREDFORLEARNER Use predictors for learners.
        %   The UsePredForLearner property is a logical matrix of size
        %   P-by-NumTrained for P predictors (columns) in the training data X and
        %   the number of trained weak learners NumTrained. An element (I,J) of
        %   this matrix is set to true if predictor I was used for training learner
        %   J and set to false otherwise.
        %
        %   For every weak learner in the ensemble, the predictors (columns) are
        %   ordered in the same way as in the training data X passed to the
        %   ensemble. For example: X has 5 columns and the 3rd column of
        %   UsePredForLearner is set to [0 1 0 1 0]'. This means that the 3rd
        %   learner in the ensemble is trained on X(:,[2 4]), not X(:,[4 2]).
        %
        %   See also classreg.learning.ensemble.CompactEnsemble.
        UsePredForLearner = [];
    end
    
    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        %NUMTRAINED Number of trained learners in the ensemble.
        %   The NumTrained property is a numeric positive scalar set to the number
        %   of trained weak learners in this ensemble.
        %
        %   See also classreg.learning.ensemble.CompactEnsemble,
        %   classreg.learning.ensemble.CompactEnsemble/Trained.
        NumTrained;
        
        %TRAINED Trained learners.
        %   The Trained property is a cell array of compact classification or
        %   regression models.
        %
        %   See also classreg.learning.ensemble.CompactEnsemble.
        Trained;
        
        %TRAINEDWEIGHTS Learner weights.
        %   The TrainedWeights is a numeric vector of weights assigned by the
        %   ensemble to its learners. The ensemble computes predicted response by
        %   aggregating weighted predictions from its learners.
        %
        %   See also classreg.learning.ensemble.CompactEnsemble,
        %   classreg.learning.ensemble.CompactEnsemble/Trained.
        TrainedWeights;
        
        %COMBINEWEIGHTS Prescription for combining weighted learner predictions.
        %   The CombineWeights property is a string describing how learner
        %   predictions are combined.
        %
        %   See also classreg.learning.ensemble.CompactEnsemble,
        %   classreg.learning.ensemble.CompactEnsemble/Trained.
        CombineWeights;
    end
    
    properties(GetAccess=public,SetAccess=protected,Hidden=true,Dependent=true)
        NTrained;
    end
    
    methods(Access=protected)
        function this = CompactEnsemble(usepredforlearner)
            this.UsePredForLearner = usepredforlearner;
        end
        
        function s = propsForDisp(this,s)
            if nargin<2 || isempty(s)
                s = struct;
            else
                if ~isstruct(s)
                    error(message('stats:classreg:learning:ensemble:CompactEnsemble:propsForDisp:BadS'));
                end
            end
            s.NumTrained = this.NumTrained;
        end
    end
    
    methods
        function t = get.NumTrained(this)
            t = length(this.Trained);
        end
        
        function t = get.NTrained(this)
            t = length(this.Trained);
        end
        
        function t = get.Trained(this)
            t = {};
            if ~isempty(this.Impl)
                t = this.Impl.Trained;
            end
        end
        
        function tw = get.TrainedWeights(this)
            tw = [];
            if ~isempty(this.Impl) && ismember('LearnerWeights',properties(this.Impl.Combiner))
                tw = this.Impl.Combiner.LearnerWeights;
            end
        end
        
        function cw = get.CombineWeights(this)
            cw = '';
            if ~isempty(this.Impl)
                cw = class(this.Impl.Combiner);
                idx = strfind(cw,'classreg.learning.combiner.');
                if ~isempty(idx)
                    cw(1:idx+length('classreg.learning.combiner.')-1) = [];
                end
            end
        end
        
        function this = set.CombineWeights(this,cw)
            if     strcmpi(cw,'weightedsum')
                this.Impl.Combiner = ...
                    classreg.learning.combiner.WeightedSum(this.TrainedWeights);
            elseif strcmpi(cw,'weightedaverage')
                this.Impl.Combiner = ...
                    classreg.learning.combiner.WeightedAverage(this.TrainedWeights);
            else
                error(message('stats:classreg:learning:ensemble:CompactEnsemble:setCombineWeights:BadArgs'));
            end
        end
    
        function this = removeLearners(this,idx)
        %REMOVELEARNERS Remove learners from ensemble.
        %   ENS=REMOVELEARNERS(ENS,IDX) returns an ensemble from which learners
        %   with indices IDX have been removed. You must pass IDX as a vector with
        %   integer values ranging from 1 to NumTrained.
        %
        %   See also classreg.learning.ensemble.CompactEnsemble, NumTrained,
        %   Trained.
        
            this.Impl = removeLearners(this.Impl,idx);
        end
    end
    
    methods(Static,Hidden)
        function score = aggregatePredict(X,combiner,trained,...
                classnames,nonzeroProbClasses,defaultScore,varargin)
            % Get sizes
            [N,D] = size(X);
            T = length(trained);

            % Any trained learners?
            % Decode input args
            [useNfort,useDfort,learnerIdx,mode] = ...
                classreg.learning.ensemble.CompactEnsemble.checkAggregateArgs(T,N,D,varargin{:});
            
            % Only the default 'ensemble' mode is allowed
            if strncmpi(mode,'ensemble',length(mode))==0
                error(message('stats:classreg:learning:ensemble:CompactEnsemble:aggregatePredict:BadMode'));
            end

            % Check useNfort
            if ~islogical(useNfort) || ~all([N T]==size(useNfort))
                error(message('stats:classreg:learning:ensemble:CompactEnsemble:aggregatePredict:UseObsForIter', N, T));
            end
            
            % Check useDfort
            if ~islogical(useDfort) || ~all([D T]==size(useDfort))
                error(message('stats:classreg:learning:ensemble:CompactEnsemble:aggregatePredict:UsePredForIter', D, T));
            end
            
            % Predict. The combiner object aggregates scores from
            % individual learners. The score returned by combiner is the
            % aggregated score.
            score = repmat(defaultScore,N,max(1,numel(classnames)));
            T = length(learnerIdx);
            if T==0
                return;
            end
            for t=1:T
                [combiner,score] = ...
                    classreg.learning.ensemble.CompactEnsemble.predictOneWithCache(...
                    combiner,X,learnerIdx(t),useNfort,useDfort,...
                    trained,classnames,nonzeroProbClasses,defaultScore);
            end
        end
            
        function v = aggregateLoss(T,X,Y,W,cost,funloss,combiner,fpredict,...
                trained,classnames,nonzeroProbClasses,scoretransform,defaultScore,varargin)
            % Decode input args
            [N,D] = size(X);
            if ~isempty(classnames)
                funloss = classreg.learning.internal.lossCheck(funloss,'classification');
            else
                funloss = classreg.learning.internal.lossCheck(funloss,'regression');
            end
            [useNfort,useDfort,learnerIdx,mode] = ...
                classreg.learning.ensemble.CompactEnsemble.checkAggregateArgs(T,N,D,varargin{:});
            
            % Any trained learners?
            if T==0
                v = NaN;
                return;
            end
            
            % If this is a classification ensemble, wrap cost inside
            % another anonymous function
            if isempty(cost)
                fwrapped = @(Y,Sfit,W) funloss(Y,Sfit,W);
            else
                fwrapped = @(Y,Sfit,W) funloss(Y,Sfit,W,cost);
            end
            
            % Get predictions for different modes. The combiner object
            % aggregates scores from individual learners. The score
            % returned by combiner is the aggregated score.
            T = length(learnerIdx);
            if     strncmpi(mode,'ensemble',length(mode))
                for t=1:T
                    [combiner,Sfit] = fpredict(combiner,X,learnerIdx(t),...
                        useNfort,useDfort,trained,classnames,nonzeroProbClasses,defaultScore);
                end
                Sfit = scoretransform(Sfit);
                v = fwrapped(Y,Sfit,W);
            elseif strncmpi(mode,'individual',length(mode))
                v = NaN(T,1);
                for t=1:T
                    [~,Sfit] = fpredict(combiner,X,learnerIdx(t),...
                        useNfort,useDfort,trained,classnames,nonzeroProbClasses,defaultScore);
                    Sfit = scoretransform(Sfit);
                    v(t) = fwrapped(Y,Sfit,W);
                end
            elseif strncmpi(mode,'cumulative',length(mode))
                v = NaN(T,1);
                for t=1:T
                    [combiner,Sfit] = fpredict(combiner,X,learnerIdx(t),...
                        useNfort,useDfort,trained,classnames,nonzeroProbClasses,defaultScore);
                    Sfit = scoretransform(Sfit);
                    v(t) = fwrapped(Y,Sfit,W);
                end
            end
        end
        
        function [combiner,score] = predictOneWithCache(combiner,X,t,...
                useNfort,useDfort,trained,classnames,nonzeroProbClasses,defaultScore)
            % Make scores of correct size
            N = size(X,1);
            K = length(classnames);
            doclass = true;
            if K==0
                doclass = false;
                K = 1;
            end
            score = NaN(N,K);
            
            % Match classes
            weak = trained{t};
            if doclass
                [~,pos] = ismember(weak.ClassSummary.ClassNames,classnames);
            end
            
            % Get scores for all ensemble members
            useNfortT = useNfort(:,t); % observations for learner t
            useDfortT = useDfort(:,t); % predictors for learner t
            if any(useNfortT)
                
                % If not all predictors have been used to train learner t,
                % treat this as proof that the ensemble has been grown by
                % random subspace. In that case, do not compute predictions
                % for observations with missing inputs. To compute
                % prediction for an observation with missing values,
                % subspace ensemble averages over learners trained on
                % non-missing inputs.
                if ~all(useDfortT)
                    nanobs = any(isnan(X(useNfortT,useDfortT)),2);
                    idxNfortT = find(useNfortT);
                    idxobs = idxNfortT(~nanobs);
                    useNfortT = false(N,1);
                    useNfortT(idxobs) = true;
                end
                    
                if doclass
                    [~,score(useNfortT,pos)] = predict(weak,X(useNfortT,useDfortT));
                else
                    score(useNfortT) = predict(weak,X(useNfortT,useDfortT));
                end
            end
            combiner = updateCache(combiner,score,t,useNfortT);
            score = cachedScore(combiner);
            
            % Assign scores only for classes with non-zero probability
            if doclass
                tf = ismember(classnames,nonzeroProbClasses);
                score(:,~tf) = defaultScore;
            end
        end
        
        function [useNfort,useDfort,learners,mode] = checkAggregateArgs(T,N,D,varargin)
            % Decode input args
            args = {'useobsforlearner' 'usepredforlearner' 'learners'     'mode'};
            defs = {                []                  []        1:T 'ensemble'};
            [useNfort,useDfort,learners,mode] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Check UseObsForIter
            if isempty(useNfort)
                useNfort = true(N,T);
            else
                if ~islogical(useNfort) || size(useNfort,1)~=N
                    error(message('stats:classreg:learning:ensemble:CompactEnsemble:checkAggregateArgs:BadUseObsForIter', N));
                end
            end
            
            % Check UsePredForIter
            if isempty(useDfort)
                useDfort = true(D,T);
            else
                if ~islogical(useDfort) || size(useDfort,1)~=D
                    error(message('stats:classreg:learning:ensemble:CompactEnsemble:checkAggregateArgs:BadUsePredForIter', D));
                end
            end            
            
            % Check learner indices and reset T
            if islogical(learners)
                if ~isvector(learners) || length(learners)~=T
                    error(message('stats:classreg:learning:ensemble:CompactEnsemble:checkAggregateArgs:BadLogicalIndices', T));
                end
                learners = find(learners);
            end
            if ~isempty(learners) && ...
                    (~isnumeric(learners) || ~isvector(learners) || min(learners)<=0 || max(learners)>T)
                error(message('stats:classreg:learning:ensemble:CompactEnsemble:checkAggregateArgs:BadNumericIndices', T));
            end
            learners = ceil(learners);
            
            % Check aggregation mode
            mode = lower(mode);
            if ~ischar(mode) || ...
                    ~any(strncmpi(mode,{'ensemble' 'individual' 'cumulative'},length(mode)))
                error(message('stats:classreg:learning:ensemble:CompactEnsemble:checkAggregateArgs:BadMode'));
            end
        end        
    end
    
end

