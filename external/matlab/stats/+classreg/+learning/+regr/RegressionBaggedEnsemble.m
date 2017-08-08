classdef RegressionBaggedEnsemble < ...
        classreg.learning.regr.RegressionEnsemble & classreg.learning.ensemble.BaggedEnsemble
%RegressionBaggedEnsemble Regression ensemble grown by resampling.
%   RegressionBaggedEnsemble combines a set of trained weak learner
%   models and data on which these learners were trained. It can predict
%   ensemble response for new data by aggregating predictions from its weak
%   learners. It also stores data used for training and can compute
%   resubstitution and out-of-bag predictions. It can resume training if
%   desired.
%
%   This class is derived from RegressionEnsemble.
%
%   RegressionBaggedEnsemble properties:
%       NumObservations       - Number of observations.
%       X                     - Matrix of predictors used to train this ensemble.
%       Y                     - Observed response used to train this ensemble.
%       W                     - Weights of observations used to train this ensemble.
%       ModelParameters       - Ensemble parameters.
%       PredictorNames        - Names of predictors used for this ensemble.
%       CategoricalPredictors - Indices of categorical predictors.
%       ResponseName          - Name of the response variable.
%       ResponseTransform     - Transformation applied to predicted regression response.
%       Method                - Ensemble algorithm used for training.
%       LearnerNames          - Names of weak learners.
%       ReasonForTermination  - Reason for stopping ensemble training.
%       FitInfo               - Ensemble fit information.
%       FitInfoDescription    - Description of ensemble fit information.
%       NumTrained            - Number of trained learners in the ensemble.
%       Trained               - Trained learners.
%       TrainedWeights        - Learner weights.
%       CombineWeights        - Prescription for combining weighted learner predictions.
%       Regularization        - Regularization results.
%       FResample             - Fraction of training data for resampling.
%       Replace               - Flag indicating if training data were sampled with replacement.
%       UseObsForLearner      - Use observations for learners.
%       UsePredForLearner     - Use predictors for learners.
%
%   RegressionBaggedEnsemble methods:
%       compact               - Compact this ensemble.
%       crossval              - Cross-validate this ensemble.
%       loss                  - Regression loss.
%       predict               - Predicted response of this model.
%       predictorImportance   - Importance of predictors for this model.
%       resubLoss             - Resubstitution regression loss.
%       resubPredict          - Resubstitution predicted response.
%       oobLoss               - Out-of-bag regression loss.
%       oobPredict            - Out-of-bag predicted response.
%       resume                - Resume training.
%       regularize            - Optimize learner weights.
%       shrink                - Discard learners with small weights.
%       cvshrink              - Cross-validate shrinking.
%
%   See also CompactRegressionEnsemble, RegressionEnsemble.

%   Copyright 2010-2013 The MathWorks, Inc.

    
    methods(Hidden)
        function this = RegressionBaggedEnsemble(X,Y,W,modelParams,...
                dataSummary,responseTransform)
            this = this@classreg.learning.regr.RegressionEnsemble(...
                X,Y,W,modelParams,dataSummary,responseTransform);
            this = this@classreg.learning.ensemble.BaggedEnsemble();
        end
    end
    
    methods(Access=protected)
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.regr.RegressionEnsemble(this,s);
            s = propsForDisp@classreg.learning.ensemble.BaggedEnsemble(this,s);
        end
    end
    
    methods(Static,Hidden)
        function this = fit(X,Y,varargin)
            error(message('stats:classreg:learning:regr:RegressionBaggedEnsemble:fit:Noop'));
        end
    end
    
    methods
        function yfit = oobPredict(this,varargin)
        %OOBPREDICT Predict out-of-bag response of the ensemble.
        %   YFIT=OOBPREDICT(ENS) returns predicted response YFIT for regression
        %   ensemble ENS for out-of-bag data. YFIT is a vector of type double with
        %   size(ENS.X,1) elements.
        %
        %   You can get indices for out-of-bag observations for weak learner L from
        %   ~ENS.UseObsForLearner(:,L).
        %
        %   YFIT=OOBPREDICT(ENS,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'learners'         - Indices of weak learners in the ensemble
        %                            ranging from 1 to NumTrained. Only these
        %                            learners are used for making predictions. By
        %                            default, all learners are used.
        %
        %   See also RegressionBaggedEnsemble, predict.        
            
            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            usenfort = ~this.ModelParams.Generator.UseObsForIter;
            yfit = predict(this,this.X,'useobsforlearner',usenfort,varargin{:});
        end
        
        function l = oobLoss(this,varargin)
        %OOBLOSS Out-of-bag regression error.
        %   L=OOBLOSS(ENS) returns mean squared error for ensemble ENS computed for
        %   out-of-bag data.
        %
        %   L=OOBLOSS(ENS,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'lossfun'   - Function handle for loss, or string representing a
        %                     built-in loss function. Available loss functions for
        %                     regression: 'mse'. If you pass a function handle FUN,
        %                     LOSS calls it as shown below:
        %                          FUN(Y,Yfit,W)
        %                     where Y, Yfit and W are numeric vectors of length N.
        %                     Y is observed response, Yfit is predicted response,
        %                     and W is observation weights. Default: 'mse'
        %       'learners'  - Indices of weak learners in the ensemble
        %                     ranging from 1 to NumTrained. Only these learners are
        %                     used for making predictions. By default, all learners
        %                     are used.
        %       'mode'      - 'ensemble' (default), 'individual' or
        %                     'cumulative'. If 'ensemble', this method returns a
        %                     scalar value for the full ensemble. If 'individual',
        %                     this method returns a vector with one element per
        %                     trained learner. If 'cumulative', this method returns
        %                     a vector in which element J is obtained by using
        %                     learners 1:J from the input list of learners.
        %
        %   See also RegressionBaggedEnsemble, loss.
            
            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            usenfort = ~this.ModelParams.Generator.UseObsForIter;
            l = loss(this,this.X,this.PrivY,'weights',this.W,...
                'useobsforlearner',usenfort,varargin{:});
        end
    end
    
end
