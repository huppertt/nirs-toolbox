classdef WeightedSum < classreg.learning.internal.DisallowVectorOps
    
%   Copyright 2010 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=protected)
        % Weights for combining scores along the 1st dimension
        LearnerWeights = []; 
        
        % Logical array of length numel(LearnerWeights) indicating which weights are
        % currently cached
        IsCached = []; 
    end
    
    properties(GetAccess=protected,SetAccess=protected)
        % Cached scores from previous call to combine()
        CachedScore = [];
    end
    
    methods
        function this = WeightedSum(learnerWeights)
            this = this@classreg.learning.internal.DisallowVectorOps();
            this.LearnerWeights = learnerWeights(:);
            this.IsCached = false(numel(this.LearnerWeights),1);
        end

        function obj = clone(this,learnerWeights)
            if nargin<2
                learnerWeights = this.LearnerWeights;
            end
            obj = classreg.learning.combiner.WeightedSum(learnerWeights);
        end

        % Empty default implementation. This method can have non-empty
        % implementations in derived classes.
        function this = addWeights(this,score,t,usenfort)
        end
        
        % Empty default implementation. This method can have non-empty
        % implementations in derived classes.
        function this = initWeights(this,score,t,usenfort)
        end
        
        function this = resetCache(this)
            this.CachedScore = [];
            this.IsCached = false(numel(this.LearnerWeights),1);
        end
        
        function score = cachedScore(this)
            score = this.CachedScore;
        end
        
        function this = updateCache(this,score,t,usenfort)
            % Get number of learners and observations
            T = length(this.LearnerWeights);
            [N K] = size(score);
            
            % Check learner index and cache size
            t = ceil(t);
            if t<1 || t>T
                error(message('stats:classreg:learning:combiner:WeightedSum:updateCache:BadLearnerIndex', T));
            end
            if ~isempty(this.CachedScore) && ~all(size(this.CachedScore)==[N K])
                error(message('stats:classreg:learning:combiner:WeightedSum:updateCache:BadCacheSize'));
            end
            
            % If first caching operation, make CachedScore.
            % Else reset combined values to zeros.
            if isempty(this.CachedScore)
                this.CachedScore = NaN(size(score));
                this.CachedScore(usenfort,:) = 0;
                this = initWeights(this,score,t,usenfort);
            else
                tf = isnan(this.CachedScore);
                if any(tf(:))
                    this.CachedScore(tf & repmat(usenfort,1,size(this.CachedScore,2))) = 0;
                end
            end

            % If the learner has been already combined or has non-positive
            % weight, do nothing
            if this.IsCached(t) || this.LearnerWeights(t)<=0
                return;
            end

            % Add
            this.CachedScore(usenfort,:) = this.CachedScore(usenfort,:) + ...
                score(usenfort,:)*this.LearnerWeights(t);
            this = addWeights(this,score,t,usenfort);
            this.IsCached(t) = true;
        end
    end
    
end
