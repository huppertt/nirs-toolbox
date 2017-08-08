classdef ClassificationEnsemble < ...
        classreg.learning.classif.FullClassificationModel & classreg.learning.ensemble.Ensemble ...
        & classreg.learning.classif.CompactClassificationEnsemble
%ClassificationEnsemble Classification ensemble.
%   ClassificationEnsemble combines a set of trained weak learner models
%   and data on which these learners were trained. It can predict ensemble
%   response for new data by aggregating predictions from its weak
%   learners. It also stores data used for training and can compute
%   resubstitution predictions. It can resume training if desired.
%
%   This class is derived from CompactClassificationEnsemble.
%
%   ClassificationEnsemble properties:
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
%       UsePredForLearner     - Use predictors for learners.
%
%   ClassificationEnsemble methods:
%       compact               - Compact this ensemble.
%       compareHoldout        - Compare two models using test data.
%       crossval              - Cross-validate this ensemble.
%       edge                  - Classification edge.
%       loss                  - Classification loss.
%       margin                - Classification margins.
%       predict               - Predicted response of this ensemble.
%       predictorImportance   - Importance of predictors for this ensemble.
%       resubEdge             - Resubstitution classification edge.
%       resubLoss             - Resubstitution classification loss.
%       resubMargin           - Resubstitution classification margins.
%       resubPredict          - Resubstitution predicted response.
%       resume                - Resume training.
%
%   See also CompactClassificationEnsemble.

%   Copyright 2010-2014 The MathWorks, Inc.

    methods(Hidden)
        function this = ClassificationEnsemble(X,Y,W,modelParams,...
                dataSummary,classSummary,scoreTransform,adjustprior)
            if nargin<8
                adjustprior = true;
            end
            if adjustprior
                [classSummary,W] = ...
                    classreg.learning.internal.adjustPrior(classSummary,Y,W);
            end
            this = this@classreg.learning.classif.FullClassificationModel(X,Y,W,modelParams,...
                dataSummary,classSummary,scoreTransform);
            this = this@classreg.learning.ensemble.Ensemble();
            this = this@classreg.learning.classif.CompactClassificationEnsemble(...
                dataSummary,classSummary,scoreTransform,[],[],[],[],[]);
            this.DefaultScore = this.ModelParams.DefaultScore;
            nlearn = this.ModelParams.NLearn/numel(this.ModelParams.LearnerTemplates);
            this = fitEnsemble(this,nlearn);
            if isa(this.ModelParams.Generator,'classreg.learning.generator.SubspaceSampler')
                this.UsePredForLearner = this.ModelParams.Generator.UsePredForIter;
            end
            
            this.DefaultScoreType = 'unknown';
            switch this.Method
                case {'AdaBoostM1' 'AdaBoostM2' 'AdaBoostMH' 'RobustBoost' ...
                        'LogitBoost' 'GentleBoost' 'RUSBoost'} % trees
                    % Non-negative for AdaBoostM2, AdaBoostMH and RUSBoost
                    this.DefaultScoreType = 'inf';
                    if ismember(this.Method,{'AdaBoostM1' 'LogitBoost' 'GentleBoost'})
                        this.TransformToProbability = ...
                            @classreg.learning.transform.doublelogit;
                        if strcmp(this.Method,'LogitBoost')
                            this.PrivContinuousLoss = @classreg.learning.loss.binodeviance;
                        else
                            this.PrivContinuousLoss = @classreg.learning.loss.exponential;
                        end
                    end
                case {'Bag' 'Subspace'}                        % discriminant and knn
                    isprob = true;
                    if this.NTrained>0
                        for t=1:this.NTrained
                            lrn = this.Trained{t};
                            if ~strcmp(lrn.ScoreType,'probability')
                                isprob = false;
                                break;
                            end
                        end
                    end
                    if isprob
                        this.DefaultScoreType = 'probability';
                        this.TransformToProbability = ...
                            @classreg.learning.transform.identity;
                        this.PrivContinuousLoss = @classreg.learning.loss.quadratic;
                    end
                otherwise                                      % 'LPBoost' 'TotalBoost'
                    this.DefaultScoreType = '01';
                    this.PrivContinuousLoss = @classreg.learning.loss.quadratic;
            end
        end
    end

    methods(Static,Hidden)
        function this = fit(X,Y,varargin)
            warning(message('stats:classreg:learning:classif:ClassificationEnsemble:fit:Noop'));
            args = {'method'};
            defs = {      ''};
            [method,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            temp = classreg.learning.FitTemplate.make(method,'type','classification',extraArgs{:});
            this = fit(temp,X,Y);
        end
    end
    
    methods
        function cmp = compact(this)
        %COMPACT Compact ensemble.
        %   CMP=COMPACT(ENS) returns an object of class
        %   CompactClassificationEnsemble holding the trained weak learners for
        %   this ensemble. The compact object does not contain X and Y used for
        %   training.
        %
        %   See also ClassificationEnsemble, CompactClassificationEnsemble.
        
            cmp = classreg.learning.classif.CompactClassificationEnsemble(...
                this.DataSummary,this.ClassSummary,...
                this.PrivScoreTransform,this.PrivScoreType,...
                this.UsePredForLearner,...
                this.DefaultScore,this.PrivContinuousLoss,this.TransformToProbability);
            if this.ModelParams.SortLearnersByWeight
                cmp.Impl = sortLearnersByWeight(this.Impl);
            else
                cmp.Impl = this.Impl;
            end
            cmp.DefaultScoreType = this.DefaultScoreType;
        end
        
        function partModel = crossval(this,varargin)
        %CROSSVAL Cross-validate this ensemble.
        %   CVENS=CROSSVAL(ENS) builds a partitioned ensemble CVENS from ensemble
        %   ENS represented by a full object for classification. You can then
        %   assess the predictive performance of this ensemble on cross-validated
        %   data using methods and properties of CVENS. By default, CVENS is built
        %   using 10-fold cross-validation on the training data. CVENS is of class
        %   ClassificationPartitionedEnsemble.
        %
        %   CVENS=CROSSVAL(ENS,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'KFold'       - Number of folds for cross-validation, a numeric
        %                       positive scalar; 10 by default.
        %       'Holdout'     - Holdout validation uses the specified
        %                       fraction of the data for test, and uses the rest of
        %                       the data for training. Specify a numeric scalar
        %                       between 0 and 1.
        %       'Leaveout'    - If 'on', use leave-one-out cross-validation.
        %       'CVPartition' - An object of class CVPARTITION; empty by default.
        %                       If a CVPARTITION object is supplied, it is used for
        %                       splitting the data into subsets.
        %       'NPrint'      - Print-out frequency, a positive integer scalar. By
        %                       default, this parameter is set to 'off' (no
        %                       print-outs). You can use this parameter to keep
        %                       track of how many cross-validation have been
        %                       trained, so far.
        %
        %   See also ClassificationEnsemble, cvpartition,
        %   classreg.learning.partition.ClassificationPartitionedEnsemble.

            partModel = crossval@classreg.learning.classif.FullClassificationModel(this,varargin{:});
        end
        
        function this = resume(this,nlearn,varargin)
        %RESUME Resume training of an ensemble.
        %   ENS=RESUME(ENS,NLEARN) trains ensemble ENS for NLEARN more cycles.
        %   NLEARN must be a positive integer scalar.
        %
        %   ENS=RESUME(ENS,NLEARN,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'NPrint' -  Print-out frequency, a positive integer scalar. By
        %                   default, this parameter is set to 'off' (no print-outs).
        %                   You can use this parameter to keep track of how many
        %                   weak learners have been trained, so far. This is useful
        %                   when you train ensembles with many learners on large
        %                   datasets.
        %
        % See also ClassificationEnsemble, fitensemble.
            
            if isempty(nlearn) || ~isnumeric(nlearn) || ~isscalar(nlearn) || nlearn<=0
                error(message('stats:classreg:learning:classif:ClassificationEnsemble:resume:BadNLearn'));
            end
            nlearn = ceil(nlearn);
            this.ModelParams.NPrint = ...
                classreg.learning.ensemble.Ensemble.checkNPrint(varargin{:});
            this.ModelParams.NLearn = this.ModelParams.NLearn + ...
                nlearn*numel(this.ModelParams.LearnerTemplates);
            this = fitEnsemble(this,nlearn);
            if isa(this.ModelParams.Generator,'classreg.learning.generator.SubspaceSampler')
                this.UsePredForLearner = this.ModelParams.Generator.UsePredForIter;
            end
        end
        
        function [labels,scores] = resubPredict(this,varargin)
        %RESUBPREDICT Predict response of the ensemble by resubstitution.
        %   [LABEL,SCORE]=RESUBPREDICT(ENS) returns class labels and scores for
        %   classification ensemble ENS for training data ENS.X. Classification
        %   labels LABEL have the same type as Y used for training. Scores SCORE
        %   are an N-by-K numeric matrix for N observations and K classes. High
        %   score value indicates that an observation likely comes from this class.
        %
        %   [LABEL,SCORE]=RESUBPREDICT(ENS,'PARAM1',val1,'PARAM2',val2,...) specifies
        %   optional parameter name/value pairs:
        %       'learners'         - Indices of weak learners in the ensemble
        %                            ranging from 1 to NumTrained. Only these
        %                            learners are used for making predictions. By
        %                            default, all learners are used.
        %
        %   See also ClassificationEnsemble, predict.
        
            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            [labels,scores] = ...
                resubPredict@classreg.learning.classif.FullClassificationModel(this,varargin{:});
        end
        
        function m = resubMargin(this,varargin)
        %RESUBMARGIN Classification margins by resubstitution.
        %   M=RESUBMARGIN(ENS) returns resubstitution classification margins.
        %   Classification margin is the difference between classification score
        %   for the true class and maximal classification score for the false
        %   classes. The returned M is a numeric column-vector of length
        %   size(ENS.X,1).
        %
        %   M=RESUBMARGIN(ENS,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'learners'         - Indices of weak learners in the ensemble
        %                            ranging from 1 to NumTrained. Only these
        %                            learners are used for making predictions. By
        %                            default, all learners are used.
        %
        %   See also ClassificationEnsemble, margin.
        
            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            m = resubMargin@classreg.learning.classif.FullClassificationModel(this,varargin{:});
        end
        
        function e = resubEdge(this,varargin)
        %RESUBEDGE Classification edge by resubstitution.
        %   E=RESUBEDGE(ENS) returns classification edge obtained by ensemble ENS
        %   for training data ENS.X and ENS.Y. Classification edge is
        %   classification margin averaged over the entire data.
        %
        %   E=RESUBEDGE(ENS,'PARAM1',val1,'PARAM2',val2,...) specifies optional
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
        %   See also ClassificationEnsemble, resubMargin, edge.
        
            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            e = resubEdge@classreg.learning.classif.FullClassificationModel(this,varargin{:});
        end
        
        function l = resubLoss(this,varargin)
        %RESUBLOSS Classification error by resubstitution.
        %   L=RESUBLOSS(ENS) returns classification error for ensemble ENS computed
        %   for training data ENS.X and ENS.Y.
        %
        %   L=RESUBLOSS(ENS,'PARAM1',val1,'PARAM2',val2,...) specifies optional
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
        %   See also ClassificationEnsemble, loss.
        
            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            l = resubLoss@classreg.learning.classif.FullClassificationModel(this,varargin{:});
        end             
    end
    
    methods(Access=protected)                
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.classif.FullClassificationModel(this,s);
            s = propsForDisp@classreg.learning.ensemble.Ensemble(this,s);
        end
        
        function this = fitEnsemble(this,nlearn)
            % Fit
            [this,trained,generator,modifier,combiner] = ...
                fitWeakLearners(this,nlearn,this.ModelParams.NPrint);            
            
            % Update generator and modifier
            this.ModelParams.Generator = generator;
            this.ModelParams.Modifier = modifier;
            
            % Make compact object
            this.Impl = classreg.learning.impl.CompactEnsembleImpl(trained,combiner);
        end
    end
    
end
