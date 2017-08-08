classdef RegressionPartitionedEnsemble < ...
        classreg.learning.partition.PartitionedEnsemble & ...
        classreg.learning.partition.RegressionPartitionedModel
%RegressionPartitionedEnsemble Cross-validated regression ensemble.
%   RegressionPartitionedEnsemble is a set of regression ensembles
%   trained on cross-validated folds. You can obtain a cross-validated
%   ensemble either by calling CROSSVAL method of a regression ensemble
%   or by calling FITENSEMBLE with one of the cross-validation options.
%
%   This class is derived from RegressionPartitionedModel.
%
%   To estimate the quality of regression by cross-validation, you can
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
%   RegressionPartitionedEnsemble properties:
%      CrossValidatedModel   - Name of the cross-validated ensemble.
%      PredictorNames        - Names of predictors used for this ensemble.
%      CategoricalPredictors - Indices of categorical predictors.
%      ResponseName          - Name of the response variable.
%      NumObservations       - Number of observations.
%      X                     - X data used to cross-validate this ensemble.
%      Y                     - Observed response values used to cross-validate this ensemble.
%      W                     - Weights of observations used to cross-validate this ensemble.
%      ModelParameters       - Cross-validation parameters.
%      Trained               - Compact regression models trained on cross-validation folds.
%      KFold                 - Number of cross-validation folds.
%      NumTrainedPerFold     - Number of trained learners per cross-validation fold.
%      Partition             - Data partition used to cross-validate this ensemble.
%      ResponseTransform     - Transformation applied to predicted regression response.
%
%   RegressionPartitionedEnsemble methods:
%      kfoldPredict          - Predict response for observations not used for training.
%      kfoldLoss             - Regression loss for observations not used for training.
%      kfoldfun              - Cross-validate function.
%      resume                - Train learners in the cross-validation folds for more cycles.
%
%   See also cvpartition, classreg.learning.regr.RegressionModel,
%   RegressionPartitionedModel.

%   Copyright 2010-2013 The MathWorks, Inc.


    methods(Hidden)
        function this = RegressionPartitionedEnsemble(X,Y,W,modelParams,...
                dataSummary,responseTransform)
            this = this@classreg.learning.partition.PartitionedEnsemble();
            this = this@classreg.learning.partition.RegressionPartitionedModel(...
                X,Y,W,modelParams,dataSummary,responseTransform);
        end
    end

    methods(Access=protected)
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.partition.PartitionedEnsemble(this,s);
            s = propsForDisp@classreg.learning.partition.RegressionPartitionedModel(this,s);
        end
    end
    
    methods(Hidden)
        function this = regularize(this,varargin)
             for t=1:numel(this.Ensemble.Trainable)
                 this.Ensemble.Trainable{t} = regularize(this.Ensemble.Trainable{t},varargin{:});
             end
         end
         
         function this = shrink(this,varargin)
             if isempty(this.Ensemble) || isempty(this.Ensemble.Impl)
                 error(message('stats:classreg:learning:partition:RegressionPartitionedEnsemble:prune:NoImpl'));
             end
             if numel(this.Ensemble.Trainable)~=numel(this.Ensemble.Trained)
                 error(message('stats:classreg:learning:partition:RegressionPartitionedEnsemble:prune:MismatchTrainedTrainableSize'));
             end
             for t=1:numel(this.Ensemble.Trainable)
                 this.Ensemble.Impl.Trained{t} = ...
                     shrink(this.Ensemble.Trainable{t},varargin{:});
             end
         end
     end
    
    methods
        function l = kfoldLoss(this,varargin)
        %KFOLDLOSS Regression loss for observations not used for training.
        %   L=KFOLDLOSS(ENS) returns loss obtained by cross-validated regression
        %   ensemble ENS. For every fold, this method computes regression loss for
        %   in-fold observations using a model trained on out-of-fold observations.
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
        %                            vector with one element per fold. if
        %                            'cumulative', this method returns a vector of
        %                            length min(NTrainedPerFold) in which element J
        %                            is obtained by averaging values across all
        %                            folds for weak learners 1:J in each fold.
        %       'lossfun'          - Function handle for loss, or string
        %                            representing a built-in loss function.
        %                            Available loss functions for regression:
        %                            'mse'. If you pass a function handle FUN, LOSS
        %                            calls it as shown below:
        %                               FUN(Y,Yfit,W)
        %                            where Y, Yfit and W are numeric vectors of
        %                            length N. Y is observed response, Yfit is
        %                            predicted response, and W is observation
        %                            weights. Default: 'mse'
        %
        %   See also RegressionPartitionedEnsemble,
        %   classreg.learning.regr.RegressionModel/loss, kfoldPredict.
            
            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            args = {                  'lossfun'};
            defs = {@classreg.learning.loss.mse};
            [funloss,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            usenfort = ~this.Ensemble.ModelParams.Generator.UseObsForIter;
            [mode,folds,partArgs] = checkEnsembleFoldArgs(this,extraArgs{:});
            if strcmp(mode,'cumulative')
                T = min(this.NTrainedPerFold(folds));
                D = numel(this.Ensemble.PredictorNames);
                trained = this.Ensemble.Impl.Trained(folds);
                useDinFold = classreg.learning.partition.PartitionedEnsemble.usePredInFold(...
                    folds,T,D,trained);
                l = classreg.learning.ensemble.CompactEnsemble.aggregateLoss(...
                    T,this.Ensemble.X,this.Ensemble.Y,this.Ensemble.W,[],funloss,this.Combiner(folds),...
                    @classreg.learning.partition.PartitionedEnsemble.predictKfoldWithCache,...
                    trained,[],[],this.Ensemble.PrivResponseTransform,NaN,...
                    'useobsforlearner',usenfort(:,folds),'mode',mode,...
                    'usepredforlearner',useDinFold,partArgs{:});
            else
                l = loss(this.Ensemble,this.Ensemble.X,this.Ensemble.Y,'lossfun',funloss,...
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
        % See also RegressionPartitionedEnsemble, fitensemble,
        % classreg.learning.regr.RegressionEnsemble/resume.
            
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
    
end
