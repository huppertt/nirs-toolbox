classdef ClassificationSVM < ...
        classreg.learning.classif.FullClassificationModel & ...
        classreg.learning.classif.CompactClassificationSVM
%ClassificationSVM Support Vector Machine model for classification.
%   ClassificationSVM is an SVM model for classification with one or two
%   classes. This model can predict response for new data. This model also
%   stores data used for training and can compute resubstitution
%   predictions.
%
%   An object of this class cannot be created by calling the constructor.
%   Use FITCSVM to create a ClassificationSVM object by fitting an SVM
%   model to training data.
%
%   This class is derived from CompactClassificationSVM.
%
%   ClassificationSVM properties:
%       NumObservations       - Number of observations.
%       X                     - Matrix of predictors used to train this model.
%       Y                     - True class labels used to train this model.
%       W                     - Weights of observations used to train this model.
%       ModelParameters       - SVM parameters.
%       PredictorNames        - Names of predictors used for this model.
%       ResponseName          - Name of the response variable.
%       ClassNames            - Names of classes in Y.
%       Cost                  - Misclassification costs.
%       Prior                 - Prior class probabilities.
%       ScoreTransform        - Transformation applied to predicted classification scores.
%       Alpha                 - Coefficients obtained by solving the dual problem.
%       Beta                  - Coefficients for the primal linear problem.
%       Bias                  - Bias term.
%       KernelParameters      - Kernel parameters.
%       Mu                    - Predictor means.
%       Sigma                 - Predictor standard deviations.
%       SupportVectors        - Support vectors.
%       SupportVectorLabels   - Support vector labels (+1 and -1).
%       BoxConstraints        - Box constraints.
%       CacheInfo             - Cache information.
%       ConvergenceInfo       - Convergence information.
%       Gradient              - Gradient values in the training data.
%       IsSupportVector       - Indices of support vectors in the training data.
%       Nu                    - Nu parameter for one-class learning.
%       NumIterations         - Number of iterations taken by optimization.
%       OutlierFraction       - Expected fraction of outliers in the training data.
%       ShrinkagePeriod       - Number of iterations between reductions of the active set.
%       Solver                - Name of the used solver.
%
%   ClassificationSVM methods:
%       compact               - Compact this model.
%       compareHoldout        - Compare two models using test data.
%       crossval              - Cross-validate this model.
%       discardSupportVectors - Discard support vectors for linear SVM.
%       edge                  - Classification edge.
%       fitPosterior          - Find transformation from SVM scores to class posterior probabilities.
%       loss                  - Classification loss.
%       margin                - Classification margins.
%       predict               - Predicted response of this model.
%       resubEdge             - Resubstitution classification edge.
%       resubLoss             - Resubstitution classification loss.
%       resubMargin           - Resubstitution classification margins.
%       resubPredict          - Resubstitution predicted response.
%       resume                - Resume training.
%
%   Example: Train an SVM model on ionosphere data.
%       load ionosphere
%       svm = fitcsvm(X,Y,'KernelFunction','gaussian','standardize',true)
%
%   See also fitcsvm, classreg.learning.classif.CompactClassificationSVM.
    
%   Copyright 2013-2014 The MathWorks, Inc.
          
    properties(SetAccess=protected,GetAccess=public,Dependent=true)
        %BOXCONSTRAINTS Box constraints.
        %   The BoxConstraints property is a vector with NumObservations positive
        %   elements. The dual coefficient alpha for observation I cannot exceed
        %   BoxConstraints(I).
        %
        %   See also fitcsvm, ClassificationSVM, Alpha, NumObservations.
        BoxConstraints;
        
        %CACHEINFO Cache information.
        %   The CacheInfo property is a struct with 2 fields:
        %       Size      - Size of memory used for caching entries of the Gram
        %                   matrix in megabytes.
        %       Algorithm - Name of the algorithm used for removing entries from
        %                   the cache when its capacity is exceeded.
        %
        %   See also fitcsvm, ClassificationSVM.
        CacheInfo;
        
        %CONVERGENCEINFO Convergence information.
        %   The ConvergenceInfo property is a struct with 10 fields:
        %       Converged                - true if optimization converged to the
        %                                  desired tolerance and false otherwise.
        %       ReasonForConvergence     - Name of the criterion used to detect
        %                                  convergence.
        %       Gap                      - Value of the attained feasibility gap
        %                                  between the dual and primal objectives.
        %       GapTolerance             - Tolerance for the feasibility gap.
        %       DeltaGradient            - Value of the attained gradient
        %                                  difference between upper and lower
        %                                  violators.
        %       DeltaGradientTolerance   - Tolerance for gradient difference
        %                                  between upper and lower violators.
        %       LargestKKTViolation      - Value of the attained largest (by
        %                                  magnitude) Karush-Kuhn-Tucker (KKT)
        %                                  violation.
        %       KKTTolerance             - Tolerance for largest KKT violation.
        %       History                  - Struct with 6 fields:
        %           NumIterations            * Array of iteration indices at which
        %                                      convergence criteria were recorded.
        %           Gap                      * Gap values at these iterations.
        %           DeltaGradient            * DeltaGradient values at these
        %                                      iterations.
        %           LargestKKTViolation      * LargestKKTViolation values at these
        %                                      iterations.
        %           NumSupportVectors        * Numbers of support vectors at these
        %                                      iterations.
        %           Objective                * Objective values at these
        %                                      iterations.
        %       Objective                - Value of the dual objective.
        %
        %   See also fitcsvm, ClassificationSVM.
        ConvergenceInfo;
        
        %GRADIENT Gradient values in the training data.
        %   The Gradient property is a vector with NumObservations elements.
        %   Element I of this vector is the value of the gradient at observation I
        %   at the end of the optimization.
        %
        %   See also fitcsvm, ClassificationSVM, NumObservations.
        Gradient;
        
        %ISSUPPORTVECTOR Indices of support vectors in the training data.
        %   The IsSupportVector property is a logical vector with NumObservations
        %   elements. IsSupportVector(I) is true if observation I is a support
        %   vector and false otherwise.
        %
        %   See also fitcsvm, ClassificationSVM, NumObservations.
        IsSupportVector;
        
        %NU Nu parameter for one-class learning.
        %   The NU property is a positive scalar. The Alpha coefficients must
        %   satisfy sum(Alpha)=NU*NumObservations.
        %
        %   See also fitcsvm, ClassificationSVM, Alpha, NumObservations.
        Nu;

        %NUMITERATIONS Number of iterations taken by optimization.
        %   The NumIterations property is an integer showing how many iterations
        %   were carried by optimization.
        %
        %   See also fitcsvm, ClassificationSVM.
        NumIterations;
        
        %OUTLIERFRACTION Expected fraction of outliers in the training data.
        %   The OutlierFraction property is a numeric scalar between 0 and 1
        %   specifying the expected fraction of outliers in the training data.
        %
        %   See also fitcsvm, ClassificationSVM.
        OutlierFraction;
        
        %SHRINKAGEPERIOD Number of iterations between reductions of the active set.
        %   The ShrinkagePeriod property is a non-negative integer specifying how
        %   often the active set was shrunk during optimization.
        %
        %   See also fitcsvm, ClassificationSVM.
        ShrinkagePeriod;
        
        %SOLVER Name of the used solver.
        %   The Solver property is a string specifying the algorithm used to solve
        %   the SVM problem.
        %
        %   See also fitcsvm, ClassificationSVM.
        Solver;
    end
    
    methods
        function a = get.BoxConstraints(this)
            a = this.Impl.C;
        end
        
        function a = get.CacheInfo(this)
            a = this.Impl.CacheInfo;
        end
        
        function a = get.ConvergenceInfo(this)
            a = this.Impl.ConvergenceInfo;
            if isfield(a,'OutlierHistory')
                a = rmfield(a,'OutlierHistory');
            end
            if isfield(a,'ChangeSetHistory')
                a = rmfield(a,'ChangeSetHistory');
            end
        end
        
        function a = get.OutlierFraction(this)
            a = this.Impl.FractionToExclude;
        end
        
        function a = get.Gradient(this)
            a = this.Impl.Gradient;
        end
        
        function a = get.Solver(this)
            a = this.ModelParams.Solver;
        end
        
        function a = get.Nu(this)
            if numel(this.ClassSummary.NonzeroProbClasses)==1
                a = this.ModelParams.Nu;
            else
                a = [];
            end
        end
        
        function a = get.NumIterations(this)
            a = this.Impl.NumIterations;
        end
                
        function a = get.IsSupportVector(this)
            a = this.Impl.IsSupportVector;
        end
        
        function a = get.ShrinkagePeriod(this)
            a = this.Impl.Shrinkage.Period;
        end
    end
    
    methods(Static,Hidden)
        function temp = template(varargin)            
            classreg.learning.FitTemplate.catchType(varargin{:});
            temp = classreg.learning.FitTemplate.make('SVM','type','classification',varargin{:});
        end
        
        function this = fit(X,Y,varargin)
            temp = ClassificationSVM.template(varargin{:});
            this = fit(temp,X,Y);
        end
    end
    
    methods(Hidden)
        function this = ClassificationSVM(X,Y,W,modelParams,...
                dataSummary,classSummary,scoreTransform)
            % Protect against being called by the user
            if nargin~=7 || ischar(W)
                error(message('stats:ClassificationSVM:ClassificationSVM:DoNotUseConstructor'));
            end
            
            % Base constructor
            this = this@classreg.learning.classif.FullClassificationModel(...
                X,Y,W,modelParams,dataSummary,classSummary,scoreTransform);
            this = this@classreg.learning.classif.CompactClassificationSVM(...
                dataSummary,classSummary,scoreTransform,[],[]);
                        
            % Remove all rows with NaNs
            nanX = any(isnan(this.X),2);
            if any(nanX)
                this.X(nanX,:)   = [];
                this.PrivY(nanX) = [];
                this.W(nanX)     = [];
            end
            if isempty(this.X)
                error(message('stats:ClassificationSVM:ClassificationSVM:NoDataAfterNaNsRemoved'));
            end
            
            % Check for alphas. The size of alphas matches the size of X. This
            % is ensured by SVMParams.
            if ~isempty(this.ModelParams.Alpha)
                this.ModelParams.Alpha(nanX) = [];
            end
            
            % If NaN values have been removed, renormalize observation weights.
            % Renormalize prior probabilities as well in case some classes now
            % have zero weights.
            if any(nanX)
                % Renormalize weights
                this.W = this.W/sum(this.W);
                
                % Save initial values. SVM always operates on the default cost
                % matrix. This is ensured by prepareData.
                nonzeroClassNames = this.ClassSummary.NonzeroProbClasses;
                prior = this.ClassSummary.Prior;
                K = numel(nonzeroClassNames);
                cost = ones(K) - eye(K);
                
                % Get matrix of class weights
                C = classreg.learning.internal.classCount(nonzeroClassNames,this.PrivY);
                WC = bsxfun(@times,C,this.W);
                Wj = sum(WC,1);
                
                % Remove observations for classes with zero prior probabilities
                % or costs
                [this.X,this.PrivY,~,WC,Wj,prior,cost,nonzeroClassNames] = ...
                    ClassificationTree.removeZeroPriorAndCost(...
                    this.X,this.PrivY,C,WC,Wj,prior,cost,nonzeroClassNames);
                                
                % Normalize priors in such a way that the priors in present
                % classes add up to one.  Normalize weights to add up to the
                % prior in the respective class.
                prior = prior/sum(prior);
                this.W = sum(bsxfun(@times,WC,prior./Wj),2);
                
                % Put processed values into summary structure
                this.ClassSummary = ...
                    classreg.learning.classif.FullClassificationModel.makeClassSummary(...
                    this.ClassSummary.ClassNames,nonzeroClassNames,prior,cost);
            end
                
            % Map classes to -1 and +1
            gidx = grp2idx(this.PrivY,this.ClassSummary.NonzeroProbClasses);
            if any(gidx==2)
                doclass = 2;
                gidx(gidx==1) = -1;
                gidx(gidx==2) = +1;
            else
                doclass = 1;
            end
                        
            this.Impl = classreg.learning.impl.SVMImpl.make(...
                this.X,gidx,this.W,...
                this.ModelParams.Alpha,...
                this.ModelParams.KernelFunction,...
                this.ModelParams.KernelPolynomialOrder,[],...
                this.ModelParams.KernelScale,this.ModelParams.KernelOffset,...
                this.ModelParams.StandardizeData,...
                doclass,...
                this.ModelParams.Solver,...
                this.ModelParams.BoxConstraint,...
                this.ModelParams.Nu,...
                this.ModelParams.IterationLimit,...
                this.ModelParams.KKTTolerance,...
                this.ModelParams.GapTolerance,...
                this.ModelParams.DeltaGradientTolerance,...
                this.ModelParams.CacheSize,...
                this.ModelParams.CachingMethod,...
                this.ModelParams.ShrinkagePeriod,...
                this.ModelParams.OutlierFraction,...
                this.ModelParams.VerbosityLevel,...
                this.ModelParams.NumPrint);
        end
    end
    
    methods(Access=protected)
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.classif.FullClassificationModel(this,s);
            s = propsForDisp@classreg.learning.classif.CompactClassificationSVM(this,s);            
            if isfield(s,'SupportVectors')
                s = rmfield(s,'SupportVectors');
            end
            if isfield(s,'SupportVectorLabels')
                s = rmfield(s,'SupportVectorLabels');
            end
            s.BoxConstraints               = this.BoxConstraints;
            s.ConvergenceInfo              = this.ConvergenceInfo;
            s.IsSupportVector              = this.IsSupportVector;
            s.Solver                       = this.Solver;
        end
    end
    
    methods
        function cmp = compact(this)
        %COMPACT Compact SVM model.
        %   CMP=COMPACT(MODEL) returns an object of class CompactClassificationSVM
        %   holding the structure of the trained SVM classifier. The compact object
        %   does not contain X and Y used for training.
        %
        %   See also fitcsvm, ClassificationSVM,
        %   classreg.learning.classif.CompactClassificationSVM.        
        
            cmp = classreg.learning.classif.CompactClassificationSVM(...
                this.DataSummary,this.ClassSummary,...
                this.PrivScoreTransform,this.PrivScoreType,...
                compact(this.Impl,this.ModelParams.SaveSupportVectors));
        end
        
        function [varargout] = resubPredict(this,varargin)
        %RESUBPREDICT Predict resubstitution response.
        %   [LABEL,SCORE]=RESUBPREDICT(OBJ) returns predicted class labels and
        %   scores for SVM model OBJ and training data OBJ.X. Classification labels
        %   LABEL have the same type as Y used for training. Scores SCORE are an
        %   N-by-K numeric matrix for N observations and K classes. The predicted
        %   label is assigned to the class with the largest score.
        %
        %   See also fitcsvm, ClassificationSVM, predict.
            
            [varargout{1:nargout}] = ...
                resubPredict@classreg.learning.classif.FullClassificationModel(this,varargin{:});
        end
        
        function [varargout] = resubLoss(this,varargin)
        %RESUBLOSS Classification error by resubstitution.
        %   L=RESUBLOSS(OBJ) returns resubstitution classification cost for model
        %   OBJ.
        %
        %   L=RESUBLOSS(DISCR,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'lossfun'   - Function handle for loss, or string representing a
        %                     built-in loss function. Available loss functions for
        %                     classification: 'binodeviance', 'classiferror',
        %                     'exponential', and 'hinge'. If you pass a function
        %                     handle FUN, LOSS calls it as shown below:
        %                          FUN(C,S,W,COST)
        %                     where C is an N-by-K logical matrix for N rows in X
        %                     and K classes in the ClassNames property, S is an
        %                     N-by-K numeric matrix, W is a numeric vector with N
        %                     elements, and COST is a K-by-K numeric matrix. C has
        %                     one true per row for the true class. S is a matrix of
        %                     posterior probabilities for classes with one row per
        %                     observation, similar to SCORE output from
        %                     PREDICT. W is a vector of observation weights. COST
        %                     is a matrix of misclassification costs. Default:
        %                     'classiferror'
        %
        %   See also fitcsvm, ClassificationSVM, ClassNames, predict, loss.

            [varargout{1:nargout}] = ...
                resubLoss@classreg.learning.classif.FullClassificationModel(this,varargin{:});
        end
        
        function this = resume(this,numIter,varargin)
        %RESUME Resume training this SVM model.
        %   MODEL=RESUME(MODEL,NUMITER) trains the SVM model MODEL for NUMITER more
        %   iterations and returns an updated model. You can resume training an SVM
        %   model if optimization has not converged and if Solver is set to 'SMO'
        %   or 'ISDA'.
        %
        %   MODEL=RESUME(MODEL,NUMITER,'PARAM1',val1,'PARAM2',val2,...) specifies
        %   optional parameter name/value pairs:
        %       'Verbose'       - Verbosity level, one of: 0, 1, or 2. Default:
        %                         Value passed to FITCSVM.
        %       'NumPrint'      - Number of iterations between consecutive
        %                         diagnostic print-outs, a non-negative integer.
        %                         RESUME uses this parameter only if you pass 1 for
        %                         the 'Verbose' parameter. Default: Value passed to
        %                         FITCSVM.
        %
        %   See also fitcsvm, ClassificationSVM, Solver.
            
            % Decode input args
            args = {                      'verbose'                'numprint'};
            defs = {this.ModelParams.VerbosityLevel this.ModelParams.NumPrint};
            [verbose,nprint,~] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            if ~isnumeric(numIter) || ~isscalar(numIter) ...
                    || isnan(numIter) || isinf(numIter) || numIter<=0
                error(message('stats:ClassificationSVM:resume:BadNumIter'));
            else
                numIter = ceil(numIter);
            end
            this.ModelParams.IterationLimit = ...
                this.ModelParams.IterationLimit + numIter;
            
            if verbose<=0
                nprint = 0;
            end
            
            gidx = grp2idx(this.PrivY,this.ClassSummary.NonzeroProbClasses);
            if any(gidx==2)
                doclass = 2;
                gidx(gidx==1) = -1;
                gidx(gidx==2) = +1;
            else
                doclass = 1;
            end
            
            this.Impl = resume(this.Impl,this.X,gidx,numIter,doclass,verbose,nprint);
        end
        
        function partModel = crossval(this,varargin)
        %CROSSVAL Cross-validate this model.
        %   CVMODEL=CROSSVAL(MODEL) builds a partitioned model CVMODEL from model
        %   MODEL represented by a full object for classification. You can then
        %   assess the predictive performance of this model on cross-validated data
        %   using methods and properties of CVMODEL. By default, CVMODEL is built
        %   using 10-fold cross-validation on the training data. CVMODEL is of
        %   class ClassificationPartitionedModel.
        %
        %   CVMODEL=CROSSVAL(MODEL,'PARAM1',val1,'PARAM2',val2,...) specifies
        %   optional parameter name/value pairs:
        %      'KFold'      - Number of folds for cross-validation, a numeric
        %                     positive scalar; 10 by default.
        %      'Holdout'    - Holdout validation uses the specified
        %                     fraction of the data for test, and uses the rest of
        %                     the data for training. Specify a numeric scalar
        %                     between 0 and 1.
        %      'Leaveout'   - If 'on', use leave-one-out cross-validation.
        %      'CVPartition' - An object of class CVPARTITION; empty by default. If
        %                      a CVPARTITION object is supplied, it is used for
        %                      splitting the data into subsets.
        %
        %   See also fitcsvm, ClassificationSVM, cvpartition,
        %   classreg.learning.partition.ClassificationPartitionedModel.

            idxBaseArg = find(ismember(lower(varargin(1:2:end)),...
                classreg.learning.FitTemplate.AllowedBaseFitObjectArgs));
            if ~isempty(idxBaseArg)
                error(message('stats:classreg:learning:classif:FullClassificationModel:crossval:NoBaseArgs', varargin{ 2*idxBaseArg - 1 }));
            end
            modelParams = this.ModelParams;
            modelParams.VerbosityLevel = 0;
            temp = classreg.learning.FitTemplate.make(this.ModelParams.Method,...
                'type','classification','scoretransform',this.PrivScoreTransform,...
                'modelparams',modelParams,'CrossVal','on',varargin{:});
            partModel = fit(temp,this.X,this.Y,'Weights',this.W,...
                'predictornames',this.DataSummary.PredictorNames,...
                'responsename',this.ResponseName,...
                'classnames',this.ClassNames,'cost',this.Cost,'prior',this.Prior);
            partModel.ScoreType = this.ScoreType;
        end
        
        function [obj,trans] = fitPosterior(obj,varargin)
        %FITPOSTERIOR Fit posterior probabilities
        %   OBJ=FITPOSTERIOR(OBJ) finds the optimal sigmoid transformation from
        %   scores to posterior probabilities for SVM model OBJ by 10-fold
        %   cross-validation. FITPOSTERIOR returns an object of type
        %   ClassificationSVM and sets the ScoreTransform property of this object
        %   to the optimal transformation. The PREDICT method for the updated model
        %   then returns posterior probabilities for 2nd output.
        %
        %   [OBJ,TRANS]=FITPOSTERIOR(OBJ) also returns TRANS, a struct with
        %   parameters of the optimal transformation from score S to posterior
        %   probability P for the positive class in OBJ.ClassNames(2). TRANS has
        %   the following fields:
        %       Type                           - String, one of: 'sigmoid', 'step'
        %                                        or 'constant'. If the two classes
        %                                        overlap, FITPOSTERIOR sets Type to
        %                                        'sigmoid'. If the two classes are
        %                                        perfectly separated, FITPOSTERIOR
        %                                        sets Type to 'step'. If one of the
        %                                        two classes has zero probability,
        %                                        FITPOSTERIOR sets Type to
        %                                        'constant.
        %          If Type is 'sigmoid', TRANS has additional fields:
        %             Slope                    - Slope A of the sigmoid
        %                                        transformation
        %                                        P(S)=1/(1+exp(A*S+B))
        %             Intercept                - Intercept B of the sigmoid
        %                                        transformation
        %                                        P(S)=1/(1+exp(A*S+B))
        %          If Type is 'step', TRANS has additional fields:
        %             PositiveClassProbability - Probability of the positive class
        %                                        in the interval between LowerBound
        %                                        and UpperBound
        %             LowerBound               - Lower bound of the interval in
        %                                        which the probability for the
        %                                        positive class is set to
        %                                        PositiveClassProbability. Below
        %                                        this bound, the probability for
        %                                        the positive class is zero.
        %             UpperBound               - Upper bound of the interval in
        %                                        which the probability for the
        %                                        positive class is set to
        %                                        PositiveClassProbability. Above
        %                                        this bound, the probability for
        %                                        the positive class is one.
        %          If Type is 'constant', TRANS has additional fields:
        %             PredictedClass           - Name of the predicted class, same
        %                                        type as OBJ.ClassNames. The
        %                                        posterior probability is one for
        %                                        this class.
        %
        %   OBJ=FITPOSTERIOR(OBJ,PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'CVPartition'           - A partition created with CVPARTITION to
        %                                 use for cross-validation.
        %       'Holdout'               - Holdout validation uses the specified
        %                                 fraction of the data for test, and uses
        %                                 the rest of the data for training.
        %                                 Specify a numeric scalar between 0 and 1.
        %       'KFold'                 - Number of folds to use for
        %                                 cross-validation, a positive integer.
        %                                 Default: 10
        %       'Leaveout'              - Use leave-one-out cross-validation by
        %                                 setting to 'on'.
        %
        %   See also fitcsvm, ClassificationSVM,
        %   classreg.learning.partition.ClassificationPartitionedModel, ClassNames,
        %   ScoreTransform.

            [obj,trans] = fitSVMPosterior(obj,varargin{:});
        end
    end
    
    methods(Static,Hidden)
        function [X,Y,W,dataSummary,classSummary,scoreTransform] = ...
                prepareData(X,Y,varargin)
            % Process input args
            args = {'classnames' 'cost' 'prior' 'scoretransform'};
            defs = {          []     []      []               []};
            [userClassNames,cost,prior,transformer,~,crArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Get class names before any rows might be removed
            allClassNames = levels(classreg.learning.internal.ClassLabel(Y));
            if isempty(allClassNames)
                error(message('stats:classreg:learning:classif:FullClassificationModel:prepareData:EmptyClassNames'));
            end
            
            % Accept real data only
            if ~isfloat(X)                
                error(message('stats:ClassificationSVM:prepareData:BadXType'));
            end
            internal.stats.checkSupportedNumeric('X',X,true);
            
            % Pre-process            
            [X,Y,W,dataSummary] = ...
                classreg.learning.FullClassificationRegressionModel.prepareDataCR(...
                X,classreg.learning.internal.ClassLabel(Y),crArgs{:});
            
            % Categorical predictors not allowed
            if ~isempty(dataSummary.CategoricalPredictors)
                error(message('stats:ClassificationSVM:prepareData:DoNotPassCategoricalPredictors'));
            end
            
            % Process class names
            [X,Y,W,userClassNames,nonzeroClassNames] = ...
                classreg.learning.classif.FullClassificationModel.processClassNames(...
                X,Y,W,userClassNames,allClassNames);
            internal.stats.checkSupportedNumeric('Weights',W,true);
          
            % Sort nonzeroClassNames to make sure -1 corresponds to 1st
            % class and +1 corresponds to 2nd class passed by the user.
            [~,loc] = ismember(userClassNames,nonzeroClassNames);
            loc(loc==0) = [];
            nonzeroClassNames = nonzeroClassNames(loc);

            % Remove missing values
            [X,Y,W] = classreg.learning.classif.FullClassificationModel.removeMissingVals(X,Y,W);
                       
            % Get matrix of class weights
            C = classreg.learning.internal.classCount(nonzeroClassNames,Y);
            WC = bsxfun(@times,C,W);
            Wj = sum(WC,1);
                      
            % Check prior
            prior = classreg.learning.classif.FullClassificationModel.processPrior(...
                prior,Wj,userClassNames,nonzeroClassNames);

            % Get costs
            cost = classreg.learning.classif.FullClassificationModel.processCost(...
                cost,prior,userClassNames,nonzeroClassNames);
        
            % Remove observations for classes with zero prior probabilities
            [X,Y,~,WC,Wj,prior,cost,nonzeroClassNames] = ...
                ClassificationTree.removeZeroPriorAndCost(...
                X,Y,C,WC,Wj,prior,cost,nonzeroClassNames);
            
            % Apply the average cost correction for 2 classes and set the
            % cost to the default value.
            if numel(nonzeroClassNames)>1
                prior = prior.*sum(cost,2)';
                cost = ones(2) - eye(2);
            end
            
            % Normalize priors in such a way that the priors in present
            % classes add up to one.  Normalize weights to add up to the
            % prior in the respective class.
            prior = prior/sum(prior);
            W = sum(bsxfun(@times,WC,prior./Wj),2);
    
            % Put processed values into summary structure
            classSummary = ...
                classreg.learning.classif.FullClassificationModel.makeClassSummary(...
                userClassNames,nonzeroClassNames,prior,cost);
                
            % Only binary and one-class SVM are allowed
            K = numel(classSummary.NonzeroProbClasses);
            if K>2
                error(message('stats:ClassificationSVM:prepareData:DoNotPassMoreThanTwoClasses'));
            end
            
            % Make output score transformation
            scoreTransform = ...
                classreg.learning.classif.FullClassificationModel.processScoreTransform(transformer);
        end        
    end
end
