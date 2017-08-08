classdef ClassificationBaggedEnsemble < ...
        classreg.learning.classif.ClassificationEnsemble & classreg.learning.ensemble.BaggedEnsemble
%ClassificationBaggedEnsemble Classification ensemble grown by resampling.
%   ClassificationBaggedEnsemble combines a set of trained weak learner
%   models and data on which these learners were trained. It can predict
%   ensemble response for new data by aggregating predictions from its weak
%   learners. It also stores data used for training and can compute
%   resubstitution and out-of-bag predictions. It can resume training if
%   desired.
%
%   This class is derived from ClassificationEnsemble.
%
%   ClassificationBaggedEnsemble properties:
%       NumObservations       - Number of observations.
%       X                     - Matrix of predictors used to train this ensemble.
%       Y                     - True class labels used to train this ensemble.
%       W                     - Weights of observations used to train this ensemble.
%       ModelParameters       - Ensemble parameters.
%       PredictorNames        - Names of predictors used for this ensemble.
%       CategoricalPredictors - Indices of categorical predictors.
%       ResponseName          - Name of the response variable.
%       ClassNames            - Names of classes in Y.
%       Cost                  - Misclassification costs.
%       Prior                 - Prior class probabilities.
%       ScoreTransform        - Transformation applied to predicted classification scores.
%       Method                - Ensemble algorithm used for training.
%       LearnerNames          - Names of weak learners.
%       ReasonForTermination  - Reason for stopping ensemble training.
%       FitInfo               - Ensemble fit information.
%       FitInfoDescription    - Description of ensemble fit information.
%       NumTrained            - Number of trained learners in the ensemble.
%       Trained               - Trained learners.
%       TrainedWeights        - Learner weights.
%       CombineWeights        - Prescription for combining weighted learner predictions.
%       FResample             - Fraction of training data for resampling.
%       Replace               - Flag indicating if training data were sampled with replacement.
%       UseObsForLearner      - Use observations for learners.
%       UsePredForLearner     - Use predictors for learners.
%
%   ClassificationBaggedEnsemble methods:
%       compact               - Compact this ensemble.
%       compareHoldout        - Compare two models using test data.
%       crossval              - Cross-validate this ensemble.
%       edge                  - Classification edge.
%       loss                  - Classification loss.
%       margin                - Classification margins.
%       predict               - Predicted response of this model.
%       predictorImportance   - Importance of predictors for this model.
%       resubEdge             - Resubstitution classification edge.
%       resubLoss             - Resubstitution classification loss.
%       resubMargin           - Resubstitution classification margins.
%       resubPredict          - Resubstitution predicted response.
%       oobEdge               - Out-of-bag classification edge.
%       oobLoss               - Out-of-bag classification loss.
%       oobMargin             - Out-of-bag classification margins.
%       oobPredict            - Out-of-bag predicted response.
%       resume                - Resume training.
%
%   See also ClassificationEnsemble, CompactClassificationEnsemble.

%   Copyright 2010-2013 The MathWorks, Inc.

    
    methods(Hidden)
        function this = ClassificationBaggedEnsemble(X,Y,W,modelParams,...
                dataSummary,classSummary,scoreTransform)
            this = this@classreg.learning.classif.ClassificationEnsemble(X,Y,W,modelParams,...
                dataSummary,classSummary,scoreTransform);
            this = this@classreg.learning.ensemble.BaggedEnsemble();
        end
    end
    
    methods(Access=protected)
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.classif.ClassificationEnsemble(this,s);
            s = propsForDisp@classreg.learning.ensemble.BaggedEnsemble(this,s);
        end
    end
    
    methods(Static,Hidden)
        function this = fit(X,Y,varargin)
            error(message('stats:classreg:learning:classif:ClassificationBaggedEnsemble:fit:Noop'));
        end
    end
        
    methods        
        function [labels,scores] = oobPredict(this,varargin)
        %OOBPREDICT Predict out-of-bag response of the ensemble.
        %   [LABEL,SCORE]=OOBPREDICT(ENS) returns class labels and scores for
        %   classification ensemble ENS for out-of-bag data. Classification labels
        %   LABEL have the same type as Y used for training. Scores SCORE are an
        %   N-by-K numeric matrix for N observations and K classes. High score
        %   value indicates that an observation likely comes from this class.
        %
        %   You can get indices for out-of-bag observations for weak learner L from
        %   ~ENS.UseObsForLearner(:,L).
        %
        %   [LABEL,SCORE]=OOBPREDICT(ENS,'PARAM1',val1,'PARAM2',val2,...) specifies
        %   optional parameter name/value pairs:
        %       'learners'         - Indices of weak learners in the ensemble
        %                            ranging from 1 to NumTrained. Only these
        %                            learners are used for making predictions. By
        %                            default, all learners are used.
        %
        %   See also ClassificationBaggedEnsemble, predict.
        
            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            usenfort = ~this.ModelParams.Generator.UseObsForIter;
            [labels,scores] = ...
                predict(this,this.X,'useobsforlearner',usenfort,varargin{:});
        end
        
        function l = oobLoss(this,varargin)
        %OOBLOSS Out-of-bag classification error.
        %   L=OOBLOSS(ENS) returns out-of-bag classification error for ensemble
        %   ENS.
        %
        %   L=OOBLOSS(ENS,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'lossfun'          - Function handle for loss, or string
        %                            representing a built-in loss function.
        %                            Available loss functions for classification:
        %                            'binodeviance', 'classiferror', and
        %                            'exponential'. If you pass a function handle
        %                            FUN, LOSS calls it as shown below:
        %                                  FUN(C,S,W,COST)
        %                            where C is an N-by-K logical matrix for N rows
        %                            in ENS.X and K classes in the ClassNames property,
        %                            S is an N-by-K numeric matrix, W is a numeric
        %                            vector with N elements, and COST is a K-by-K
        %                            numeric matrix. C has one true per row for the
        %                            true class. S is a matrix of predicted scores
        %                            for classes with one row per observation,
        %                            similar to SCORE output from PREDICT. W is a
        %                            vector of observation weights. COST is a
        %                            matrix of misclassification costs. Default:
        %                            'classiferror'
        %       'learners'         - Indices of weak learners in the ensemble
        %                            ranging from 1 to NumTrained. Only these
        %                            learners are used for making predictions. By
        %                            default, all learners are used.
        %       'mode'             - 'ensemble' (default), 'individual' or
        %                            'cumulative'. If 'ensemble', this method
        %                            returns a scalar value for the full ensemble.
        %                            If 'individual', this method returns a vector
        %                            with one element per trained learner. If
        %                            'cumulative', this method returns a vector in
        %                            which element J is obtained by using learners
        %                            1:J from the input list of learners.
        %
        %   See also ClassificationBaggedEnsemble, loss.
            
            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            usenfort = ~this.ModelParams.Generator.UseObsForIter;
            l = loss(this,this.X,this.PrivY,'weights',this.W,...
                'useobsforlearner',usenfort,varargin{:});
        end
        
        function m = oobMargin(this,varargin)
        %OOBMARGIN Out-of-bag classification margins.
        %   M=OOBMARGIN(ENS) returns out-of-bag classification margins.
        %   Classification margin is the difference between classification score
        %   for the true class and maximal classification score for the false
        %   classes. The returned M is a numeric column-vector of length
        %   size(ENS.X,1).
        %
        %   You can get indices for out-of-bag observations for weak learner L from
        %   ~ENS.UseObsForLearner(:,L).
        %
        %   M=OOBMARGIN(ENS,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'learners'         - Indices of weak learners in the ensemble
        %                            ranging from 1 to NumTrained. Only these
        %                            learners are used for making predictions. By
        %                            default, all learners are used.
        %
        %   See also ClassificationBaggedEnsemble, margin.
            
            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            usenfort = ~this.ModelParams.Generator.UseObsForIter;
            m = margin(this,this.X,this.PrivY,...
                'useobsforlearner',usenfort,varargin{:});
        end
        
        function e = oobEdge(this,varargin)
        %OOBEDGE Out-of-bag classification edge.
        %   E=OOBEDGE(ENS) returns out-of-bag classification edge for ensemble ENS.
        %   Classification edge is classification margin averaged over the entire
        %   data.
        %
        %   You can get indices for out-of-bag observations for weak learner L from
        %   ~ENS.UseObsForLearner(:,L).
        %
        %   E=OOBEDGE(ENS,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'learners'         - Indices of weak learners in the ensemble
        %                            ranging from 1 to NumTrained. Only these
        %                            learners are used for making predictions. By
        %                            default, all learners are used.
        %       'mode'             - 'ensemble' (default), 'individual' or
        %                            'cumulative'. If 'ensemble', this method
        %                            returns a scalar value for the full ensemble.
        %                            If 'individual', this method returns a vector
        %                            with one element per trained learner. If
        %                            'cumulative', this method returns a vector in
        %                            which element J is obtained by using learners
        %                            1:J from the input list of learners.
        %
        %   See also ClassificationBaggedEnsemble, oobMargin, edge.
            
            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            usenfort = ~this.ModelParams.Generator.UseObsForIter;
            e = edge(this,this.X,this.PrivY,'weights',this.W,...
                'useobsforlearner',usenfort,varargin{:});
        end
    end
    
end

