classdef ClassificationPartitionedModel < classreg.learning.partition.PartitionedModel
%ClassificationPartitionedModel Cross-validated classification model.
%   ClassificationPartitionedModel is a set of classification models
%   trained on cross-validated folds. You can obtain a cross-validated
%   model either by calling CROSSVAL method of a classification model or by
%   training a classification model with one of the cross-validation
%   options.
%
%   To estimate the quality of classification by cross-validation, you can
%   use KFOLD methods. Every KFOLD method uses models trained on in-fold
%   observations to predict response for out-of-fold observations. For
%   example, you cross-validate using 5 folds. In this case, every training
%   fold contains roughly 4/5 of data and every test fold contains roughly
%   1/5 of data. The first model stored in Trained{1} was trained on X and
%   Y with the first 1/5 excluded, the second model stored in Trained{2}
%   was trained on X and Y with the second 1/5 excluded and so on. When you
%   call KFOLDPREDICT, it computes predictions for the first 1/5 of the
%   data using the first model, for the second 1/5 of data using the second
%   model and so on. In short, response for every observation is computed by
%   KFOLDPREDICT using the model trained without this observation.
%
%   ClassificationPartitionedModel properties:
%      CrossValidatedModel   - Name of the cross-validated model.
%      PredictorNames        - Names of predictors used for this model.
%      CategoricalPredictors - Indices of categorical predictors.
%      ResponseName          - Name of the response variable.
%      NumObservations       - Number of observations.
%      X                     - X data used to cross-validate this model.
%      Y                     - True class labels used to cross-validate this model.
%      W                     - Weights of observations used to cross-validate this model.
%      ModelParameters       - Cross-validation parameters.
%      Trained               - Compact classifiers trained on cross-validation folds.
%      KFold                 - Number of cross-validation folds.
%      Partition             - Data partition used to cross-validate this model.
%      ClassNames            - Names of classes in Y.
%      Cost                  - Misclassification costs.
%      Prior                 - Prior class probabilities.
%      ScoreTransform        - Transformation applied to predicted classification scores.
%
%   ClassificationPartitionedModel methods:
%      kfoldPredict          - Predict response for observations not used for training.
%      kfoldLoss             - Classification loss for observations not used for training.
%      kfoldMargin           - Classification margins for observations not used for training.
%      kfoldEdge             - Classification edge for observations not used for training.
%      kfoldfun              - Cross-validate function.
%
%   See also cvpartition, classreg.learning.classif.ClassificationModel.

%   Copyright 2010-2014 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        %CLASSNAMES Names of classes in Y.
        %   The ClassNames property is an array containing the class names for the
        %   response variable Y.
        %
        %   See also ClassificationPartitionedModel.
        ClassNames;
    end
        
    properties(GetAccess=public,SetAccess=public,Dependent=true)
        %COST Misclassification costs.
        %   The Cost property is a square matrix with misclassification costs.
        %   Cost(I,J) is the cost of misclassifying class ClassNames(I) as class
        %   ClassNames(J).
        %
        %   See also ClassificationPartitionedModel.
        Cost;
        
        %PRIOR Prior class probabilities.
        %   The Prior property is a vector with prior probabilities for classes.
        %
        %   See also ClassificationPartitionedModel.
        Prior;
        
        %SCORETRANSFORM Transformation applied to predicted classification scores.
        %   The ScoreTransform property is a string describing how raw
        %   classification scores predicted by the model are transformed. You can
        %   assign a function handle or one of the following strings to this
        %   property: 'none', 'doublelogit', 'identity', 'invlogit', 'ismax',
        %   'logit', 'sign', 'symmetricismax', 'symmetriclogit', and 'symmetric'.
        %   You can use either 'identity' or 'none' for the identity
        %   transformation.
        %
        %   See also ClassificationPartitionedModel.
        ScoreTransform;
    end
       
    properties(GetAccess=public,SetAccess=public,Dependent=true,Hidden=true)
        %See ClassificationModel for help.
        ScoreType;
    end
    
    properties(GetAccess=public,SetAccess=protected,Dependent=true,Hidden=true)
        %See ClassificationModel for help.
        ContinuousLoss;
    end
    
    properties(GetAccess=public,SetAccess=public,Hidden=true)
        DefaultLoss = @classreg.learning.loss.classiferror;
        LabelPredictor = @classreg.learning.classif.ClassificationModel.maxScore;
    end
    
    properties(GetAccess=public,SetAccess=protected,Hidden=true)
        PrivScoreTransform = [];
        PrivScoreType = [];
    end
    
    methods
        function cnames = get.ClassNames(this)
            cnames = this.Ensemble.ClassNames;
        end
        
        function cost = get.Cost(this)
            cost = this.Ensemble.Cost;
        end
        
        function this = set.Cost(this,cost)
            this.Ensemble = setLearnersCost(this.Ensemble,cost);
        end
        
        function prior = get.Prior(this)
            prior = this.Ensemble.Prior;
        end
                
        function this = set.Prior(this,prior)
            this.Ensemble = setLearnersPrior(this.Ensemble,prior);
        end
        
        function st = get.ScoreTransform(this)
            st = classreg.learning.internal.convertScoreTransform(...
                this.PrivScoreTransform,'string',[]);
        end
        
        function this = set.ScoreTransform(this,st)
            this.PrivScoreTransform = ...
                classreg.learning.internal.convertScoreTransform(st,...
                'handle',numel(this.Ensemble.ClassSummary.ClassNames));
            this.PrivScoreType = [];
        end
        
        function st = get.ScoreType(this)
            st = getScoreType(this);
        end
        
        function this = set.ScoreType(this,st)
            this.PrivScoreType = classreg.learning.internal.convertScoreType(st);
        end
        
        function cl = get.ContinuousLoss(this)
            cl = getContinuousLoss(this);
        end
    end
    
    methods(Hidden)
        function this = ClassificationPartitionedModel(X,Y,W,modelParams,...
                dataSummary,classSummary,scoreTransform)
            this = this@classreg.learning.partition.PartitionedModel();
            this.Ensemble = classreg.learning.classif.ClassificationEnsemble(...
                X,Y,W,modelParams,dataSummary,classSummary,scoreTransform,false);
            this.PrivScoreTransform = this.Ensemble.PrivScoreTransform;
            if this.Ensemble.NTrained>0
                this.DefaultLoss = this.Ensemble.Trained{1}.DefaultLoss;
                this.LabelPredictor = this.Ensemble.Trained{1}.LabelPredictor;
            end
        end
    end
        
    methods(Access=protected)
        function st = getScoreType(this)
            if this.Ensemble.NTrained>0 ...
                    && isequal(this.PrivScoreTransform,@classreg.learning.transform.identity)
                st = this.Ensemble.Trained{1}.ScoreType;
            elseif ~isempty(this.PrivScoreType)
                st = this.PrivScoreType;
            else
                st = 'unknown';
            end
        end
        
        function cl = getContinuousLoss(this)
            cl = [];
            if this.Ensemble.NTrained>0
                if     isequal(this.PrivScoreTransform,@classreg.learning.transform.identity)
                    cl = this.Ensemble.Trained{1}.ContinuousLoss;
                elseif strcmp(this.ScoreType,'probability')
                    cl = @classreg.learning.loss.quadratic;
                end
            end            
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.partition.PartitionedModel(this,s);
            cnames = this.ClassNames;
            if ischar(cnames)
                s.ClassNames = cnames;
            else
                s.ClassNames = cnames';
            end
            s.ScoreTransform = this.ScoreTransform;
        end
    end
    
    methods
        function [varargout] = kfoldPredict(this,varargin)
        %KFOLDPREDICT Predict response for observations not used for training.
        %   [LABEL,SCORE]=KFOLDPREDICT(OBJ) returns class labels and scores
        %   predicted by cross-validated classification model OBJ. Classification
        %   labels LABEL have the same type as Y used for training. Scores SCORE
        %   are an N-by-K numeric matrix for N observations and K classes. For
        %   every fold, this method predicts class labels and scores for in-fold
        %   observations using a model trained on out-of-fold observations.
        %
        %   [LABEL,SCORE,COST]=KFOLDPREDICT(OBJ) also returns an N-by-K matrix of
        %   misclassification costs COST for some cross-validated models. If
        %   PREDICT method of a model returns misclassification costs, KFOLDPREDICT
        %   method of the cross-validated model returns misclassification costs as
        %   well.
        %
        %   See also ClassificationPartitionedModel,
        %   classreg.learning.classif.ClassificationModel/predict.

            [~,score] = kfoldPredict@classreg.learning.partition.PartitionedModel(this,varargin{:});
            [varargout{1:nargout}] = this.LabelPredictor(this.ClassNames,...
                this.Prior,this.Cost,score,this.PrivScoreTransform);
        end
        
        function m = kfoldMargin(this,varargin)
        %KFOLDMARGIN Classification margins for observations not used for training.
        %   M=KFOLDMARGIN(OBJ) returns classification margins obtained by
        %   cross-validated classification model OBJ. For every fold, this method
        %   computes classification margins for in-fold observations using a model
        %   trained on out-of-fold observations.
        %
        %   See also ClassificationPartitionedModel,
        %   classreg.learning.classif.ClassificationModel/margin, kfoldPredict.
            
            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            classreg.learning.partition.PartitionedModel.catchFolds(varargin{:});
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            [mode,~,args] = checkFoldArgs(this,varargin{:});
            usenfort = ~this.Ensemble.ModelParams.Generator.UseObsForIter;
            m = margin(this.Ensemble,this.Ensemble.X,this.Ensemble.PrivY,...
                'useobsforlearner',usenfort,'mode',mode,args{:});
        end
        
        function e = kfoldEdge(this,varargin)
        %KFOLDEDGE Classification edge for observations not used for training.
        %   E=KFOLDEDGE(OBJ) returns classification edge (average classification
        %   margin) obtained by cross-validated classification model OBJ. For every
        %   fold, this method computes classification edge for in-fold observations
        %   using a model trained on out-of-fold observations.
        %
        %   L=KFOLDEDGE(OBJ,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'Folds' - Indices of folds ranging from 1 to KFold. Use only these
        %                 folds for predictions. By default, all folds are used.
        %       'Mode'  - 'average' (default) or 'individual'. If 'average', this
        %                 method averages over all folds. If 'individual', this
        %                 method returns a vector with one element per fold.
        %
        %   See also ClassificationPartitionedModel,
        %   classreg.learning.classif.ClassificationModel/edge, kfoldMargin.

            e = kfoldLoss(this,'lossfun',@classreg.learning.loss.classifedge,varargin{:});
        end

        function l = kfoldLoss(this,varargin)
        %KFOLDLOSS Classification loss for observations not used for training.
        %   L=KFOLDLOSS(OBJ) returns loss obtained by cross-validated
        %   classification model OBJ. For every fold, this method computes
        %   classification loss for in-fold observations using a model trained on
        %   out-of-fold observations.
        %
        %   L=KFOLDLOSS(OBJ,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'Folds'            - Indices of folds ranging from 1 to KFold. Use
        %                            only these folds for predictions. By default,
        %                            all folds are used.
        %       'Mode'             - 'average' (default) or 'individual'. If
        %                            'average', this method averages over all
        %                             folds. If 'individual', this method returns a
        %                             vector with one element per fold.
        %       'LossFun'          - Function handle for loss, or string
        %                            representing a built-in loss function.
        %                            Available loss functions for classification:
        %                            'binodeviance', 'classiferror', 'exponential',
        %                            and 'mincost'. If you pass a function handle
        %                            FUN, LOSS calls it as shown below:
        %                                  FUN(C,S,W,COST)
        %                            where C is an N-by-K logical matrix for N rows
        %                            in X and K classes in the ClassNames property,
        %                            S is an N-by-K numeric matrix, W is a numeric
        %                            vector with N elements, and COST is a K-by-K
        %                            numeric matrix. C has one true per row for the
        %                            true class. S is a matrix of predicted scores
        %                            for classes with one row per observation,
        %                            similar to SCORE output from PREDICT. W is a
        %                            vector of observation weights. COST is a
        %                            matrix of misclassification costs. Default:
        %                            Same as for the class of OBJ.Trained{1}. See
        %                            help for LOSS method of that class.
        %
        %   See also ClassificationPartitionedModel,
        %   classreg.learning.classif.ClassificationModel/loss, kfoldPredict.
            
            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            [mode,folds,extraArgs] = checkFoldArgs(this,varargin{:});

            usenfort = ~this.Ensemble.ModelParams.Generator.UseObsForIter;

            X = this.Ensemble.X;
            Y = this.Ensemble.PrivY;
            W = this.Ensemble.W;
        
            args = {       'lossfun'};
            defs = {this.DefaultLoss};
            [funloss,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,extraArgs{:});
            
            [X,C,W,~,usenfort] = prepareDataForLoss(this.Ensemble,X,Y,W,usenfort);
            l = classreg.learning.ensemble.CompactEnsemble.aggregateLoss(...
                this.Ensemble.NTrained,X,C,W,this.Ensemble.Cost,funloss,...
                this.Ensemble.Impl.Combiner,@classreg.learning.ensemble.CompactEnsemble.predictOneWithCache,...
                this.Ensemble.Impl.Trained,this.Ensemble.ClassSummary.ClassNames,this.Ensemble.ClassSummary.NonzeroProbClasses,...
                this.PrivScoreTransform,this.Ensemble.DefaultScore,...
                'useobsforlearner',usenfort,'mode',mode,'learners',folds,extraArgs{:});
        end
    end
    
    methods(Static,Hidden)
        function this = loadobj(obj)
            if isempty(obj.PrivScoreTransform)
                this = classreg.learning.partition.ClassificationPartitionedModel(...
                    obj.Ensemble.X,obj.Ensemble.PrivY,obj.Ensemble.W,...
                    obj.Ensemble.ModelParams,...
                    obj.Ensemble.DataSummary,obj.Ensemble.ClassSummary,...
                    obj.Ensemble.PrivScoreTransform);
            else
                % Load 11b or later
                this = obj;
            end
        end
    end
    
end
