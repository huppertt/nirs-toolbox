classdef WeightedAverage < classreg.learning.combiner.WeightedSum

%   Copyright 2010 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=protected)
        CachedWeights = [];
    end
    
    methods
        function this = WeightedAverage(learnerWeights)
            this = this@classreg.learning.combiner.WeightedSum(learnerWeights);
        end

        function obj = clone(this,learnerWeights)
            if nargin<2
                learnerWeights = this.LearnerWeights;
            end
            obj = classreg.learning.combiner.WeightedAverage(learnerWeights);
        end
        
        function score = cachedScore(this)
            score = bsxfun(@rdivide,this.CachedScore,this.CachedWeights);
        end
        
        function this = addWeights(this,score,t,usenfort)
            this.CachedWeights(usenfort) = this.CachedWeights(usenfort) + ...
                this.LearnerWeights(t);
        end
        
        function this = initWeights(this,score,t,usenfort)
            this.CachedWeights = zeros(size(score,1),1);
        end
    end

end
