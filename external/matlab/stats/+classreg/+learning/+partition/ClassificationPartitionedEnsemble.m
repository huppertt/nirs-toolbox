classdef ClassificationPartitionedEnsemble < ...
        classreg.learning.partition.PartitionedEnsemble & ...
        classreg.learning.partition.ClassificationPartitionedModel
%ClassificationPartitionedEnsemble Cross-validated classification ensemble.
%   ClassificationPartitionedEnsemble is a set of classification ensembles
%   trained on cross-validated folds. You can obtain a cross-validated
%   ensemble either by calling CROSSVAL method of a classification ensemble
%   or by calling FITENSEMBLE with one of the cross-validation options.
%
%   This class is derived from ClassificationPartitionedModel.
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
%   ClassificationPartitionedEnsemble properties:
%      CrossValidatedModel   - Name of the cross-validated ensemble.
%      PredictorNames        - Names of predictors used for this ensemble.
%      CategoricalPredictors - Indices of categorical predictors.
%      ResponseName          - Name of the response variable.
%      NumObservations       - Number of observations.
%      X                     - X data used to cross-validate this ensemble.
%      Y                     - True class labels used to cross-validate this ensemble.
%      W                     - Weights of observations used to cross-validate this ensemble.
%      ModelParameters       - Cross-validation parameters.
%      Trained               - Compact classifiers trained on cross-validation folds.
%      KFold                 - Number of cross-validation folds.
%      NumTrainedPerFold     - Number of trained learners per cross-validation fold.
%      Partition             - Data partition used to cross-validate this ensemble.
%      ClassNames            - Names of classes in Y.
%      Cost                  - Misclassification costs.
%      Prior                 - Prior class probabilities.
%      ScoreTransform        - Transformation applied to predicted classification scores.
%
%   ClassificationPartitionedEnsemble methods:
%      kfoldPredict          - Predict response for observations not used for training.
%      kfoldLoss             - Classification loss for observations not used for training.
%      kfoldMargin           - Classification margins for observations not used for training.
%      kfoldEdge             - Classification edge for observations not used for training.
%      kfoldfun              - Cross-validate function.
%      resume                - Train learners in the cross-validation folds for more cycles.
%
%   See also cvpartition, classreg.learning.classif.ClassificationModel,
%   ClassificationPartitionedModel.

%   Copyright 2010-2014 The MathWorks, Inc.

    
    methods(Hidden)
        function this = ClassificationPartitionedEnsemble(X,Y,W,modelParams,...
                dataSummary,classSummary,scoreTransform)
            this = this@classreg.learning.partition.PartitionedEnsemble();
            this = this@classreg.learning.partition.ClassificationPartitionedModel(...
                X,Y,W,modelParams,dataSummary,classSummary,scoreTransform);
        end
    end
       
    methods(Access=protected)        
        function scoreType = getScoreType(this)
            scoreType = getScoreType@classreg.learning.partition.ClassificationPartitionedModel(this);
            if this.Ensemble.NTrained>0 ...
                    && isequal(this.PrivScoreTransform,this.Ensemble.Trained{1}.TransformToProbability)
                scoreType = 'probability';
            end
        end
        
        function cl = getContinuousLoss(this)
            cl = [];
            if this.Ensemble.NTrained>0
                if     isequal(this.PrivScoreTransform,@classreg.learning.transform.identity)
                    cl = this.Ensemble.Trained{1}.ContinuousLoss;
                elseif isequal(this.PrivScoreTransform,this.Ensemble.Trained{1}.TransformToProbability)
                    cl = @classreg.learning.loss.quadratic;
                end
            end            
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.partition.PartitionedEnsemble(this,s);
            s = propsForDisp@classreg.learning.partition.ClassificationPartitionedModel(this,s);
        end
    end
    
    methods
        function e = kfoldEdge(this,varargin)
        %KFOLDEDGE Classification edge for observations not used for training.
        %   E=KFOLDEDGE(ENS) returns classification edge (average classification
        %   margin) obtained by cross-validated classification ensemble ENS. For every
        %   fold, this method computes classification edge for in-fold observations
        %   using an ensemble trained on out-of-fold observations.
        %
        %   L=KFOLDEDGE(ENS,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'folds' - Indices of folds ranging from 1 to KFold. Use only these
        %                 folds for predictions. By default, all folds are used.
        %       'mode'  - 'average' (default), 'individual' or 'cumulative'. If
        %                 'average', this method returns a scalar value averaged
        %                 over all folds. If 'individual', this method returns a
        %                 vector of length KFold with one element per fold. if
        %                 'cumulative', this method returns a vector of length
        %                 min(NumTrainedPerFold) in which element J is obtained by
        %                 averaging values across all folds for weak learners 1:J
        %                 in each fold.
        %
        %   See also ClassificationPartitionedEnsemble,
        %   classreg.learning.classif.ClassificationModel/edge, kfoldMargin.

            e = kfoldLoss(this,'lossfun',@classreg.learning.loss.classifedge,varargin{:});
        end
        
        function l = kfoldLoss(this,varargin)
        %KFOLDLOSS Classification loss for observations not used for training.
        %   L=KFOLDLOSS(ENS) returns loss obtained by cross-validated
        %   classification ensemble ENS. For every fold, this method computes
        %   classification loss for in-fold observations using a model trained on
        %   out-of-fold observations.
        %
        %   L=KFOLDLOSS(ENS,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'folds'            - Indices of folds ranging from 1 to KFold. Use
        %                            only these folds for predictions. By default,
        %                            all folds are used.
        %       'mode'             - 'average' (default), 'individual' or
        %                            'cumulative'. If 'average', this method
        %                            returns a scalar value averaged over all
        %                            folds. If 'individual', this method returns a
        %                            vector of length KFold with one element per
        %                            fold. if 'cumulative', this method returns a
        %                            vector of length min(NumTrainedPerFold) in
        %                            which element J is obtained by averaging
        %                            values across all folds for weak learners 1:J
        %                            in each fold.
        %       'lossfun'          - Function handle for loss, or string
        %                            representing a built-in loss function.
        %                            Available loss functions for classification:
        %                            'binodeviance', 'classiferror', and
        %                            'exponential'. If you pass a function handle
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
        %                            'classiferror'
        %
        %   See also ClassificationPartitionedEnsemble,
        %   classreg.learning.classif.ClassificationModel/loss, kfoldPredict.

            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            args = {                           'lossfun'};
            defs = {@classreg.learning.loss.classiferror};
            [funloss,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            C = classreg.learning.internal.classCount(...
                this.Ensemble.ClassSummary.ClassNames,this.Ensemble.PrivY);
            usenfort = ~this.Ensemble.ModelParams.Generator.UseObsForIter;
            [mode,folds,partArgs] = checkEnsembleFoldArgs(this,extraArgs{:});
            if strcmp(mode,'cumulative')
                T = min(this.NTrainedPerFold(folds));
                D = numel(this.Ensemble.PredictorNames);
                trained = this.Ensemble.Impl.Trained(folds);
                useDinFold = classreg.learning.partition.PartitionedEnsemble.usePredInFold(...
                    folds,T,D,trained);
                l = classreg.learning.ensemble.CompactEnsemble.aggregateLoss(...
                    T,this.Ensemble.X,C,this.Ensemble.W,this.Ensemble.Cost,...
                    funloss,this.Combiner(folds),...
                    @classreg.learning.partition.PartitionedEnsemble.predictKfoldWithCache,...
                    trained,this.Ensemble.ClassSummary.ClassNames,...
                    this.Ensemble.ClassSummary.NonzeroProbClasses,this.Ensemble.PrivScoreTransform,...
                    this.Ensemble.DefaultScore,...
                    'useobsforlearner',usenfort(:,folds),'mode',mode,...
                    'usepredforlearner',useDinFold,partArgs{:});
            else
                l = loss(this.Ensemble,this.Ensemble.X,this.Ensemble.PrivY,'lossfun',funloss,...
                    'weights',this.Ensemble.W,'useobsforlearner',usenfort,'mode',mode,...
                    'learners',folds,partArgs{:});
            end
        end
        
        function this = resume(this,nlearn,varargin)
        %RESUME Resume training learners on cross-validation folds.
        %   ENS=RESUME(ENS,NLEARN) trains an ensemble in every fold for NLEARN more
        %   cycles and returns an updated cross-validated ensemble ENS. You must
        %   pass NLEARN as a positive integer scalar.
        %
        %   ENS=RESUME(ENS,NLEARN,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'NPrint' -  Print-out frequency, a positive integer scalar. By
        %                   default, this parameter is set to 'off' (no print-outs).
        %                   You can use this parameter to keep track of how many
        %                   folds have been trained with additional cycles, so far.
        %                   For example, if you set 'NPrint' to 1, RESUME will
        %                   display a message after every fold is updated.
        %
        % See also ClassificationPartitionedEnsemble, fitensemble,
        % classreg.learning.classif.ClassificationEnsemble/resume.

            nprint = classreg.learning.ensemble.Ensemble.checkNPrint(varargin{:});
            trainable = resumePartitionedWithPrint(this,nlearn,nprint);
            T = numel(trainable);
            trained = cell(T,1);
            for t=1:T
                trained{t} = compact(trainable{t});
            end
            this.Ensemble.Trainable = trainable;
            this.Ensemble.Impl.Trained = trained;
        end
    end
    
    methods(Static,Hidden)
        function this = loadobj(obj)
            if isempty(obj.PrivScoreTransform)
                this = classreg.learning.partition.ClassificationPartitionedEnsemble(...
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
