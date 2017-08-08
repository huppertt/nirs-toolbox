classdef ClassificationPartitionedECOC < classreg.learning.partition.ClassificationPartitionedModel
%ClassificationPartitionedECOC Cross-validated error-correcting output code (ECOC) model.
%   ClassificationPartitionedECOC is a set of ECOC models trained on
%   cross-validated folds. You can obtain a cross-validated model either by
%   calling CROSSVAL method of a ClassificationECOC model or by training a
%   classification model with one of the cross-validation options.
%
%   To estimate the quality of classification by cross-validation, you can
%   use KFOLD methods. Every KFOLD method uses models trained on in-fold
%   observations to predict responses for out-of-fold observations. For
%   example, you cross-validate using 5 folds. In this case, every training
%   fold contains roughly 4/5 of data and every test fold contains roughly
%   1/5 of data. The first model stored in Trained{1} was trained on X and
%   Y with the first 1/5 excluded, the second model stored in Trained{2}
%   was trained on X and Y with the second 1/5 excluded and so on. When you
%   call KFOLDPREDICT, it computes predictions for the first 1/5 of the
%   data using the first model, for the second 1/5 of data using the second
%   model and so on. In short, a response for every observation is computed
%   by KFOLDPREDICT using the model trained without this observation.
%
%   ClassificationPartitionedECOC properties:
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
%      BinaryY               - Class labels for binary learners.
%      BinaryLoss            - Default binary loss function for prediction.
%      CodingMatrix          - Coding matrix.
%
%   ClassificationPartitionedECOC methods:
%      kfoldPredict          - Predict response for observations not used for training.
%      kfoldLoss             - Classification loss for observations not used for training.
%      kfoldMargin           - Classification margins for observations not used for training.
%      kfoldEdge             - Classification edge for observations not used for training.
%      kfoldfun              - Cross-validate function.
%
%   See also cvpartition, fitcecoc, ClassificationECOC.

%   Copyright 2014 The MathWorks, Inc.

    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        %BINARYY Class labels for binary learners.
        %   If the same coding matrix is used across all folds, the BinaryY
        %   property is an N-by-L matrix for N observations and L binary learners
        %   specified by the columns of the coding matrix. Its elements take values
        %   -1, 0 or +1. If element (I,J) of this matrix is -1, observation I is
        %   included in the negative class for binary learner J; if this element is
        %   +1, observation I is included in the positive class for binary learner
        %   J; and if this element is 0, observation I is not used for training
        %   binary learner J.
        %
        %   If the coding matrix varies across the folds, the BinaryY property is
        %   empty.
        %
        %   See also classreg.learning.partition.ClassificationPartitionedECOC,
        %   CodingMatrix.
        BinaryY;
        
        %BINARYLOSS Default binary loss function for prediction.
        %   The BinaryLoss property is a string specifying the default function for
        %   computing loss incurred by each binary learner.
        %
        %   See also classreg.learning.partition.ClassificationPartitionedECOC,
        %   kfoldPredict.
        BinaryLoss;
        
        %CODINGMATRIX Coding matrix.
        %   If the same coding matrix is used across all folds, the CodingMatrix
        %   property is a K-by-L matrix for K classes and L binary learners. Its
        %   elements take values -1, 0 or +1. If element (I,J) of this matrix is
        %   -1, class I is included in the negative class for binary learner J; if
        %   this element is +1, class I is included in the positive class for
        %   binary learner J; and if this element is 0, class I is not used for
        %   training binary learner J.
        %
        %   If the coding matrix varies across the folds, the CodingMatrix property
        %   is empty. In this case, use the Trained property to get the coding
        %   matrix for each fold. For example, OBJ.Trained{1}.CodingMatrix returns
        %   the coding matrix in the first fold of the cross-validated ECOC model
        %   OBJ.
        %
        %   See also classreg.learning.partition.ClassificationPartitionedECOC,
        %   BinaryY.
        CodingMatrix;
    end
    
    methods
        function bY = get.BinaryY(this)
            M = this.CodingMatrix;            
            L = size(M,2);
            N = this.NObservations;
            bY = zeros(N,L);
            
            for l=1:L
                neg = M(:,l)==-1;
                pos = M(:,l)==1;
                isneg = ismember(this.Ensemble.PrivY,...
                    this.Ensemble.ClassSummary.ClassNames(neg));
                ispos = ismember(this.Ensemble.PrivY,...
                    this.Ensemble.ClassSummary.ClassNames(pos));
                bY(isneg,l) = -1;
                bY(ispos,l) =  1;
            end
        end
        
        function bl = get.BinaryLoss(this)
            learners = this.Ensemble.Trained;
            if isempty(learners)
                bl = '';
                return;
            end
            bl = learners{1}.BinaryLoss;
        end
        
        function M = get.CodingMatrix(this)
            M = [];
            learners = this.Ensemble.Trained;
            T = numel(learners);
            
            if T==0
                return;
            end
            
            M1 = [];
            
            for t=1:T
                if ~isempty(learners{t})
                    M = learners{t}.CodingMatrix;
                    [~,pos] = ismember(this.Ensemble.ClassSummary.ClassNames,...
                        learners{t}.ClassSummary.ClassNames);
                    M = M(pos,:);
                    
                    if isempty(M1)
                        M1 = M;
                    else
                        if ~isequal(M1,M)
                            M = [];
                            return;
                        end
                    end
                end
            end            
        end
    end
    
    methods(Hidden)
        function this = ClassificationPartitionedECOC(X,Y,W,modelParams,...
                dataSummary,classSummary,scoreTransform)
            this = this@classreg.learning.partition.ClassificationPartitionedModel(...
                X,Y,W,modelParams,dataSummary,classSummary,scoreTransform);
            this.DefaultLoss = @classreg.learning.loss.classiferror;
            this.LabelPredictor = @classreg.learning.classif.ClassificationModel.maxScore;
        end        
    end
    
    methods
        function [labels,negloss,pscore,posterior] = kfoldPredict(this,varargin)
        %KFOLDPREDICT Predict response for observations not used for training.
        %   [LABEL,NEGLOSS]=KFOLDPREDICT(OBJ) returns class labels LABEL and
        %   negated values of the average binary loss per class NEGLOSS predicted
        %   by cross-validated ECOC model OBJ. Classification labels LABEL have the
        %   same type as Y used for training. Negative loss values NEGLOSS are an
        %   N-by-K matrix for N observations (rows) in OBJ.X and K classes in
        %   OBJ.ClassNames. For every fold, this method predicts class labels and
        %   negative loss values for in-fold observations using a model trained on
        %   out-of-fold observations. The predicted label is assigned to the class
        %   with the largest negated average binary loss, or equivalently smallest
        %   average binary loss.
        %
        %   [~,~,PBSCORE]=KFOLDPREDICT(OBJ) also returns an N-by-L matrix of
        %   positive-class scores predicted by the binary learners for N
        %   observations in OBJ.X and L binary learners in OBJ.BinaryLearners. If
        %   the coding matrix varies across the folds, KFOLDPREDICT returns PBSCORE
        %   as an empty array.
        %
        %   [~,~,~,POSTERIOR]=KFOLDPREDICT(OBJ) also returns posterior probability
        %   estimates, an N-by-K matrix for N observations in OBJ.X and K classes
        %   in OBJ.ClassNames. KFOLDPREDICT cannot compute these estimates unless
        %   you passed 'FitPosterior' as true to FITCECOC. If you set
        %   'FitPosterior' to false for FITCECOC and request 4th output,
        %   KFOLDPREDICT throws an error.
        %
        %   [...]=KFOLDPREDICT(OBJ,X,'PARAM1',val1,'PARAM2',val2,...) specifies
        %   optional parameter name/value pairs:
        %       'BinaryLoss'           - Function handle, or string representing a
        %                                built-in function for computing loss
        %                                induced by each binary learner. Available
        %                                loss functions for binary learners with
        %                                scores in the (-inf,+inf) range are:
        %                                'hamming', 'linear', 'exponential',
        %                                'binodeviance' and 'hinge'. Available loss
        %                                functions for binary learners with scores
        %                                in the [0,1] range are: 'hamming' and
        %                                'quadratic'. If you pass a function handle
        %                                FUN, KFOLDPREDICT calls it as shown below:
        %                                      FUN(M,F)
        %                                where M is a K-by-L coding matrix saved in
        %                                the CodingMatrix property and F is a
        %                                1-by-L row-vector of scores computed by
        %                                the binary learners. Default: Value of the
        %                                BinaryLoss property
        %       'Decoding'             - String specifying the decoding scheme, either
        %                                'lossbased' or 'lossweighted'. Default:
        %                                'lossweighted'
        %       'NumKLInitializations' - Non-negative integer specifying the number
        %                                of random initial guesses for fitting
        %                                posterior probabilities by minimization of
        %                                the Kullback-Leibler divergence.
        %                                KFOLDPREDICT ignores this parameter unless
        %                                you request 4th output and set
        %                                'PosteriorMethod' to 'kl'. Default: 0
        %       'Options'              - A struct that contains options specifying
        %                                whether to use parallel computation. This
        %                                argument can be created by a call to
        %                                STATSET. Set 'Options' to
        %                                statset('UseParallel',true) to use
        %                                parallel computation.
        %       'PosteriorMethod'      - String specifying how posterior
        %                                probabilities are fitted, either 'qp' or
        %                                'kl'. If 'qp', multiclass probabilities
        %                                are fitted by solving a least-squares
        %                                problem by quadratic programming. The 'qp'
        %                                method requires an Optimization Toolbox
        %                                license. If 'kl', multiclass probabilities
        %                                are fitted by minimizing the
        %                                Kullback-Leibler divergence between the
        %                                predicted and expected posterior
        %                                probabilities returned by the binary
        %                                learners. KFOLDPREDICT ignores this
        %                                parameter unless you request 4th output.
        %                                Default: 'kl'
        %       'Verbose'              - Non-negative integer specifying the
        %                                verbosity level, either 0 or 1.
        %                                KFOLDPREDICT does not display any
        %                                diagnostic messages at verbosity level 0
        %                                and displays diagnostic messages at
        %                                verbosity level 1. Default: 0
        %
        %   See also classreg.learning.partition.ClassificationPartitionedECOC,
        %   classreg.learning.classif.CompactClassificationECOC/predict, fitcecoc,
        %   BinaryLoss, CodingMatrix, statset.
            
            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            % Do not attempt to catch the 'folds' argument because
            % kfoldLoss (which needs this argument) directly calls
            % kfoldPredict.
            % classreg.learning.partition.PartitionedModel.catchFolds(varargin{:});
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            [~,folds,args] = checkFoldArgs(this,varargin{:});
            
            N = size(this.X,1);
            K = numel(this.Ensemble.ClassSummary.ClassNames);            
            negloss = NaN(N,K);

            pscore = [];
            if nargout>2
                M = this.CodingMatrix;
                if isempty(M)
                    warning(message('stats:classreg:learning:partition:ClassificationPartitionedECOC:kfoldPredict:CodingMatrixSizeVaries'));
                else
                    L = size(M,2);
                    pscore = NaN(N,L);
                end
            end
            
            posterior = [];
            if nargout>3
                posterior = NaN(N,K);
            end
            
            learners = this.Ensemble.Trained;
            uofl = ~this.Ensemble.ModelParams.Generator.UseObsForIter;
            
            for k=1:numel(folds)
                t = folds(k);
                if ~isempty(learners{t})
                    iuse = uofl(:,t);
                    if     nargout<3 || (nargout<4 && isempty(pscore))
                        [~,negloss(iuse,:)] = ...
                            predict(learners{t},this.X(iuse,:),args{:});
                    elseif nargout<4
                        [~,negloss(iuse,:),pscore(iuse,:)] = ...
                            predict(learners{t},this.X(iuse,:),args{:}); %#ok<AGROW>
                    else
                        if isempty(pscore)
                            [~,negloss(iuse,:),~,posterior(iuse,:)] = ...
                                predict(learners{t},this.X(iuse,:),args{:}); %#ok<AGROW>
                        else
                            [~,negloss(iuse,:),pscore(iuse,:),posterior(iuse,:)] = ...
                                predict(learners{t},this.X(iuse,:),args{:}); %#ok<AGROW>
                        end
                    end
                end
            end
            
            labels = this.LabelPredictor(this.ClassNames,...
                this.Prior,this.Cost,negloss,@(x)x);
        end
                
        function l = kfoldLoss(this,varargin)
        %KFOLDLOSS Classification error for observations not used for training.
        %   L=KFOLDLOSS(OBJ) returns error obtained by cross-validated
        %   classification model OBJ. For every fold, this method computes
        %   classification error for in-fold observations using a model trained on
        %   out-of-fold observations.
        %
        %   L=KFOLDLOSS(OBJ,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'Folds'                - Indices of folds ranging from 1 to KFold.
        %                                Use only these folds for predictions. By
        %                                default, all folds are used.
        %       'Mode'                 - 'average' (default) or 'individual'. If
        %                                'average', this method averages over all
        %                                folds. If 'individual', this method
        %                                returns a vector with one element per
        %                                fold.
        %       'LossFun'              - Function handle for loss, or string
        %                                'classiferror' representing a built-in
        %                                loss function. If you pass a function
        %                                handle FUN, KFOLDLOSS calls it as shown
        %                                below:
        %                                         FUN(C,S,W,COST)
        %                                where C is an N-by-K logical matrix for N
        %                                rows in X and K classes in the ClassNames
        %                                property, S is an N-by-K numeric matrix, W
        %                                is a numeric vector with N elements, and
        %                                COST is a K-by-K numeric matrix. C has one
        %                                true per row for the true class. S is a
        %                                matrix of negated loss values for classes
        %                                with one row per observation, similar to
        %                                NEGLOSS output from KFOLDPREDICT. W is a
        %                                vector of observation weights. COST is a
        %                                matrix of misclassification costs.
        %                                Default: 'classiferror'
        %       'BinaryLoss'           - Function handle, or string representing a
        %                                built-in function for computing loss
        %                                induced by each binary learner. Available
        %                                loss functions for binary learners with
        %                                scores in the (-inf,+inf) range are:
        %                                'hamming', 'linear', 'exponential',
        %                                'binodeviance' and 'hinge'. Available loss
        %                                functions for binary learners with scores
        %                                in the [0,1] range are: 'hamming' and
        %                                'quadratic'. If you pass a function handle
        %                                FUN, KFOLDPREDICT calls it as shown below:
        %                                      FUN(M,F)
        %                                where M is a K-by-L coding matrix saved in
        %                                the CodingMatrix property and F is a
        %                                1-by-L row-vector of scores computed by
        %                                the binary learners. Default: Value of the
        %                                BinaryLoss property
        %       'Decoding'             - String specifying the decoding scheme, either
        %                                'lossbased' or 'lossweighted'. Default:
        %                                'lossweighted'
        %       'Options'              - A struct that contains options specifying
        %                                whether to use parallel computation. This
        %                                argument can be created by a call to
        %                                STATSET. Set 'Options' to
        %                                statset('UseParallel',true) to use
        %                                parallel computation.
        %       'Verbose'              - Non-negative integer specifying the
        %                                verbosity level, either 0 or 1.
        %                                KFOLDLOSS does not display any
        %                                diagnostic messages at verbosity level 0
        %                                and displays diagnostic messages at
        %                                verbosity level 1. Default: 0
        %
        %   See also classreg.learning.partition.ClassificationPartitionedECOC,
        %   kfoldPredict.
            
        % The call to kfoldLoss@ClassificationPartitionedModel won't do
        % because CompactEnsemble knows nothing about the 'binaryloss'
        % argument.
%             l = kfoldLoss@classreg.learning.partition.ClassificationPartitionedModel(...
%                 this,varargin{:});
        
            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            [mode,folds,extraArgs] = checkFoldArgs(this,varargin{:});

            % Process input args
            args = {       'lossfun'};
            defs = {this.DefaultLoss};
            [funloss,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,extraArgs{:});
            
            % Check loss function
            funloss = classreg.learning.internal.lossCheck(funloss,'classification');

            % Get negative loss values
            [~,negloss] = kfoldPredict(this,'folds',folds,extraArgs{:});
            
            % This is just to apply ScoreTransform to negloss.
            [~,negloss] = this.LabelPredictor(...
                this.ClassNames,this.Prior,this.Cost,negloss,@(x)x);

            % Class counts
            C = classreg.learning.internal.classCount(...
                this.Ensemble.ClassSummary.ClassNames,this.Ensemble.PrivY);

            % Apply loss function to predictions
            W = this.W;
            if     strncmpi(mode,'ensemble',length(mode))
                l = funloss(C,negloss,W,this.Cost);
            elseif strncmpi(mode,'individual',length(mode))
                uofl = ~this.Ensemble.ModelParams.Generator.UseObsForIter;
                T = numel(folds);
                l = NaN(T,1);
                for k=1:T
                    t = folds(k);
                    iuse = uofl(:,t);
                    l(k) = funloss(C(iuse,:),negloss(iuse,:),W(iuse),this.Cost);
                end
            else
                error(message('stats:classreg:learning:partition:PartitionedModel:checkFoldArgs:BadMode'));
            end
        end
    
        function m = kfoldMargin(this,varargin)
        %KFOLDMARGIN Classification margins for observations not used for training.
        %   M=KFOLDMARGIN(OBJ) returns classification margins obtained by
        %   cross-validated classification model OBJ. For every fold, this method
        %   computes classification margins for in-fold observations using a model
        %   trained on out-of-fold observations.
        %
        %   M=KFOLDMARGIN(OBJ,X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies
        %   optional parameter name/value pairs:
        %       'BinaryLoss'           - Function handle, or string representing a
        %                                built-in function for computing loss
        %                                induced by each binary learner. Available
        %                                loss functions for binary learners with
        %                                scores in the (-inf,+inf) range are:
        %                                'hamming', 'linear', 'exponential',
        %                                'binodeviance' and 'hinge'. Available loss
        %                                functions for binary learners with scores
        %                                in the [0,1] range are: 'hamming' and
        %                                'quadratic'. If you pass a function handle
        %                                FUN, KFOLDPREDICT calls it as shown below:
        %                                      FUN(M,F)
        %                                where M is a K-by-L coding matrix saved in
        %                                the CodingMatrix property and F is a
        %                                1-by-L row-vector of scores computed by
        %                                the binary learners. Default: Value of the
        %                                BinaryLoss property
        %       'Decoding'             - String specifying the decoding scheme, either
        %                                'lossbased' or 'lossweighted'. Default:
        %                                'lossweighted'
        %       'Options'              - A struct that contains options specifying
        %                                whether to use parallel computation. This
        %                                argument can be created by a call to
        %                                STATSET. Set 'Options' to
        %                                statset('UseParallel',true) to use
        %                                parallel computation.
        %       'Verbose'              - Non-negative integer specifying the
        %                                verbosity level, either 0 or 1.
        %                                KFOLDMARGIN does not display any
        %                                diagnostic messages at verbosity level 0
        %                                and displays diagnostic messages at
        %                                verbosity level 1. Default: 0
        %
        %   See also classreg.learning.partition.ClassificationPartitionedECOC,
        %   ClassificationECOC/margin, kfoldPredict.
        
            classreg.learning.ensemble.Ensemble.catchUOFL(varargin{:});
            classreg.learning.partition.PartitionedModel.catchFolds(varargin{:});
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            
            [~,negloss] = kfoldPredict(this,varargin{:});
            
            C = classreg.learning.internal.classCount(...
                this.Ensemble.ClassSummary.ClassNames,this.Ensemble.PrivY);
            m = classreg.learning.loss.classifmargin(C,negloss);
        end
        
        function e = kfoldEdge(this,varargin)
        %KFOLDEDGE Classification edge for observations not used for training.
        %   E=KFOLDEDGE(OBJ) returns classification edge (average classification
        %   margin) obtained by cross-validated classification model OBJ. For every
        %   fold, this method computes classification edge for in-fold observations
        %   using a model trained on out-of-fold observations.
        %
        %   E=KFOLDEDGE(OBJ,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'Folds'               - Indices of folds ranging from 1 to KFold.
        %                               Use only these folds for predictions. By
        %                               default, all folds are used.
        %       'Mode'                - 'average' (default) or 'individual'. If
        %                               'average', this method averages over all
        %                               folds. If 'individual', this method returns
        %                               a vector with one element per fold.
        %       'BinaryLoss'           - Function handle, or string representing a
        %                                built-in function for computing loss
        %                                induced by each binary learner. Available
        %                                loss functions for binary learners with
        %                                scores in the (-inf,+inf) range are:
        %                                'hamming', 'linear', 'exponential',
        %                                'binodeviance' and 'hinge'. Available loss
        %                                functions for binary learners with scores
        %                                in the [0,1] range are: 'hamming' and
        %                                'quadratic'. If you pass a function handle
        %                                FUN, KFOLDPREDICT calls it as shown below:
        %                                         FUN(M,F)
        %                                where M is a K-by-L coding matrix saved in
        %                                the CodingMatrix property and F is a
        %                                1-by-L row-vector of scores computed by
        %                                the binary learners. Default: Value of the
        %                                BinaryLoss property
        %       'Decoding'             - String specifying the decoding scheme, either
        %                                'lossbased' or 'lossweighted'. Default:
        %                                'lossweighted'
        %       'Options'              - A struct that contains options specifying
        %                                whether to use parallel computation. This
        %                                argument can be created by a call to
        %                                STATSET. Set 'Options' to
        %                                statset('UseParallel',true) to use
        %                                parallel computation.
        %       'Verbose'              - Non-negative integer specifying the
        %                                verbosity level, either 0 or 1. KFOLDEDGE
        %                                does not display any diagnostic messages
        %                                at verbosity level 0 and displays
        %                                diagnostic messages at verbosity level 1.
        %                                Default: 0
        %
        %   See also classreg.learning.partition.ClassificationPartitionedModel,
        %   ClassificationECOC/edge, kfoldMargin.            
            
            e = kfoldEdge@classreg.learning.partition.ClassificationPartitionedModel(this,varargin{:});
        end
    end
    
end
