classdef ClassificationECOC < ...
        classreg.learning.classif.FullClassificationModel & ...
        classreg.learning.classif.CompactClassificationECOC
%ClassificationECOC Error-correcting output code (ECOC) model.
%   ClassificationECOC is an ECOC model for multiclass learning. This model
%   can predict response for new data. This model also stores data used for
%   training and can compute resubstitution predictions.
%
%   An object of this class cannot be created by calling the constructor.
%   Use FITCECOC to create a ClassificationECOC object by fitting an ECOC
%   model to training data.
%
%   This class is derived from CompactClassificationECOC.
%
%   ClassificationECOC properties:
%       NumObservations       - Number of observations.
%       X                     - Matrix of predictors used to train this model.
%       Y                     - True class labels used to train this model.
%       W                     - Weights of observations used to train this model.
%       ModelParameters       - ECOC parameters.
%       PredictorNames        - Names of predictors used for this model.
%       CategoricalPredictors - Indices of categorical predictors.
%       ResponseName          - Name of the response variable.
%       ClassNames            - Names of classes in Y.
%       Cost                  - Misclassification costs.
%       Prior                 - Prior class probabilities.
%       ScoreTransform        - Transformation applied to predicted classification scores.
%       BinaryLearners        - Binary learners.
%       BinaryLoss            - Default binary loss function for prediction.
%       CodingMatrix          - Coding matrix.
%       LearnerWeights        - Weights for binary learners.
%       BinaryY               - Class labels for binary learners.
%       CodingName            - Name of the coding design.
%
%   ClassificationECOC methods:
%       compact               - Compact this model.
%       compareHoldout        - Compare two models using test data.
%       crossval              - Cross-validate this model.
%       discardSupportVectors - Discard support vectors for linear SVM.
%       edge                  - Classification edge.
%       loss                  - Classification loss.
%       margin                - Classification margins.
%       predict               - Predicted response of this model.
%       resubEdge             - Resubstitution classification edge.
%       resubLoss             - Resubstitution classification loss.
%       resubMargin           - Resubstitution classification margins.
%       resubPredict          - Resubstitution predicted response.
%
%   Example: Train a "one versus one" SVM model on Fisher iris data
%       load fisheriris
%       ecoc = fitcecoc(meas,species)
%
%   See also fitcecoc, classreg.learning.classif.CompactClassificationECOC.
    
%   Copyright 2014 The MathWorks, Inc.

    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        %BINARYY Class labels for binary learners.
        %   The BinaryY property is an N-by-L matrix for N observations and L
        %   binary learners specified by the columns of the coding matrix. Its
        %   elements take values -1, 0 or +1. If element (I,J) of this matrix is
        %   -1, observation I is included in the negative class for binary learner
        %   J; if this element is +1, observation I is included in the positive
        %   class for binary learner J; and if this element is 0, observation I is
        %   not used for training binary learner J.
        %
        %   See also ClassificationECOC, CodingMatrix.
        BinaryY;
        
        %CODINGNAME Name of the coding design.
        %   The CodingName property is a string, one of: 'onevsone', 'onevsall',
        %   'binarycomplete', 'ternarycomplete', 'ordinal', 'sparserandom',
        %   'denserandom' or 'custom'.
        %
        %   See also ClassificationECOC, CodingMatrix.
        CodingName;
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
                isneg = ismember(this.PrivY,this.ClassSummary.ClassNames(neg));
                ispos = ismember(this.PrivY,this.ClassSummary.ClassNames(pos));
                bY(isneg,l) = -1;
                bY(ispos,l) =  1;
            end
        end
        
        function dn = get.CodingName(this)
            if ischar(this.ModelParams.Coding)
                dn = this.ModelParams.Coding;
            else
                dn = 'custom';
            end
        end
    end
    
    methods(Static,Hidden)
        function temp = template(varargin)
            classreg.learning.FitTemplate.catchType(varargin{:});
            temp = classreg.learning.FitTemplate.make('ECOC','type','classification',varargin{:});
        end
        
        function this = fit(X,Y,varargin)
            temp = ClassificationECOC.template(varargin{:});
            this = fit(temp,X,Y);
        end
    end
    
    methods(Hidden)
        function this = ClassificationECOC(X,Y,W,modelParams,...
                dataSummary,classSummary,scoreTransform)
            % Protect against being called by the user
            if nargin~=7 || ischar(W)
                error(message('stats:ClassificationECOC:ClassificationECOC:DoNotUseConstructor'));
            end
                        
            % Base constructor
            this = this@classreg.learning.classif.FullClassificationModel(...
                X,Y,W,modelParams,dataSummary,classSummary,scoreTransform);
            this = this@classreg.learning.classif.CompactClassificationECOC(...
                dataSummary,classSummary,scoreTransform,[],[],[]);
            
            % Make a coding matrix
            K = numel(this.ClassSummary.ClassNames);
            if ischar(this.ModelParams.Coding)
                M = designecoc(K,this.ModelParams.Coding);
            else
                M = this.ModelParams.Coding;
            end
                        
            % Check learners. This check should be kept in the constructor
            % because the number of columns in the coding matrix can vary
            % for the two random designs. If the number of columns was
            % fixed, this check could be imposed in
            % ECOCParams/fillDefaultParams once for all cross-validation
            % folds.
            L = size(M,2);
            learners = this.ModelParams.BinaryLearners;
            if iscell(learners)
                if numel(learners)~=L
                    error(message('stats:ClassificationECOC:ClassificationECOC:BadNumberOfLearners',...
                        numel(learners),L));
                end
            else
                learners = repmat({learners},L,1);
            end
            
            C = classreg.learning.internal.classCount(...
                this.ClassSummary.ClassNames,this.PrivY);
            
            [this.BinaryLearners,this.LearnerWeights] = ...
                localFitECOC(learners,M,this.X,C,this.W,...
                this.Prior,this.Cost,...
                this.ModelParams.FitPosterior,...
                this.ModelParams.Options,...
                this.ModelParams.VerbosityLevel);
            
            this.CodingMatrix = M;
            
            [this.BinaryLoss,this.DefaultScoreType] = ...
                classreg.learning.classif.CompactClassificationECOC.analyzeLearners(...
                this.BinaryLearners);
        end
        
    end
      
    methods(Access=protected)
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.classif.CompactClassificationECOC(this,s);
            if isfield(s,'CodingMatrix')
                s = rmfield(s,'CodingMatrix');
            end
            s.CodingName = this.CodingName;
         end
    end
    
    methods        
        function partModel = crossval(this,varargin)
        %CROSSVAL Cross-validate this model.
        %   CVMODEL=CROSSVAL(MODEL) builds a partitioned model CVMODEL from model
        %   MODEL represented by a full object for classification. You can then
        %   assess the predictive performance of this model on cross-validated data
        %   using methods and properties of CVMODEL. By default, CVMODEL is built
        %   using 10-fold cross-validation on the training data. CVMODEL is of
        %   class ClassificationPartitionedECOC.
        %
        %   CVMODEL=CROSSVAL(MODEL,'PARAM1',val1,'PARAM2',val2,...) specifies
        %   optional parameter name/value pairs:
        %      'KFold'       - Number of folds for cross-validation, a numeric
        %                      positive scalar; 10 by default.
        %      'Holdout'     - Holdout validation uses the specified
        %                      fraction of the data for test, and uses the rest of
        %                      the data for training. Specify a numeric scalar
        %                      between 0 and 1.
        %      'Leaveout'    - If 'on', use leave-one-out cross-validation.
        %      'CVPartition' - An object of class CVPARTITION; empty by default. If
        %                      a CVPARTITION object is supplied, it is used for
        %                      splitting the data into subsets.
        %      'Options'     - A struct that contains options specifying
        %                      whether to use parallel computation. This argument
        %                      can be created by a call to STATSET. Set 'Options'
        %                      to statset('UseParallel',true) to use parallel
        %                      computation. 
        %
        %   See also fitcsvm, ClassificationECOC, cvpartition, statset,
        %   classreg.learning.partition.ClassificationPartitionedECOC.

            % Catch base arguments
            idxBaseArg = find(ismember(lower(varargin(1:2:end)),...
                classreg.learning.FitTemplate.AllowedBaseFitObjectArgs));
            if ~isempty(idxBaseArg)
                error(message('stats:classreg:learning:classif:FullClassificationModel:crossval:NoBaseArgs', varargin{ 2*idxBaseArg - 1 }));
            end
            
            % Catch useful arguments
            args = {'options'};
            defs = {       []};
            [paropts,~,extraArgs] = internal.stats.parseArgs(args,defs,varargin{:});
            if ~isempty(paropts) && ~isstruct(paropts)
                error(message('stats:ClassificationECOC:crossval:BadOptionsType'));
            end
            
            % Get model parameters and make a template. 
            % Verbosity level is always minimal.
            modelParams = this.ModelParams;
            modelParams.VerbosityLevel = 0;
            temp = classreg.learning.FitTemplate.make(this.ModelParams.Method,...
                'type','classification','modelparams',modelParams,...
                'crossval','on','options',paropts,extraArgs{:});
            
            % Fit from template
            partModel = fit(temp,this.X,this.Y,'Weights',this.W,...
                'predictornames',this.DataSummary.PredictorNames,...
                'categoricalpredictors',this.CategoricalPredictors,...
                'responsename',this.ResponseName,...
                'classnames',this.ClassNames,'cost',this.Cost,'prior',this.Prior);
        end
        
        function cmp = compact(this)            
        %COMPACT Compact ECOC model.
        %   CMP=COMPACT(MODEL) returns an object of class CompactClassificationECOC
        %   holding the structure of the trained ECOC classifier. The compact
        %   object does not contain X and Y used for training.
        %
        %   See also fitcecoc, ClassificationECOC,
        %   classreg.learning.classif.CompactClassificationECOC.
            
            cmp = classreg.learning.classif.CompactClassificationECOC(...
                this.DataSummary,this.ClassSummary,this.PrivScoreTransform,...
                this.BinaryLearners,this.LearnerWeights,this.CodingMatrix);
        end
                
        function [varargout] = resubPredict(this,varargin)
        %RESUBPREDICT Predict response of the ECOC model.
        %   [LABEL,NEGLOSS]=RESUBPREDICT(MODEL) returns predicted class labels
        %   LABEL and negated values of the average binary loss per class NEGLOSS
        %   for ECOC model MODEL and training data MODEL.X. Classification labels
        %   LABEL have the same type as Y used for training. Negative loss values
        %   NEGLOSS are an N-by-K matrix for N observations (rows) in MODEL.X and K
        %   classes in MODEL.ClassNames. The predicted label is assigned to the
        %   class with the largest negated average binary loss, or equivalently
        %   smallest average binary loss.
        %
        %   [~,~,PBSCORE] = RESUBPREDICT(MODEL) also returns positive-class
        %   scores predicted by the binary learners, an N-by-L matrix for N
        %   observations in MODEL.X and L binary learners in MODEL.BinaryLearners.
        %
        %   [~,~,~,POSTERIOR]=RESUBPREDICT(MODEL) also returns posterior
        %   probability estimates, an N-by-K matrix for N observations in MODEL.X
        %   and K classes in MODEL.ClassNames. RESUBPREDICT cannot compute these
        %   estimates unless you passed 'FitPosterior' as true to FITCECOC. If you
        %   set 'FitPosterior' to false for FITCECOC and request 4th output,
        %   RESUBPREDICT throws an error.
        %
        %   [...]=RESUBPREDICT(MODEL,'PARAM1',val1,'PARAM2',val2,...) specifies
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
        %                                FUN, PREDICT calls it as shown below:
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
        %                                RESUBPREDICT ignores this parameter unless
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
        %                                learners. RESUBPREDICT ignores this
        %                                parameter unless you request 4th output.
        %                                Default: 'kl'
        %       'Verbose'              - Non-negative integer specifying the
        %                                verbosity level, either 0 or 1.
        %                                RESUBPREDICT does not display any
        %                                diagnostic messages at verbosity level 0
        %                                and displays diagnostic messages at
        %                                verbosity level 1. Default: 0
        %
        %   See also ClassificationECOC, predict.
        
            [varargout{1:nargout}] = ...
                resubPredict@classreg.learning.classif.FullClassificationModel(this,varargin{:});
        end
        
        function [varargout] = resubLoss(this,varargin)
        %RESUBLOSS Classification error by resubstitution.
        %   L=RESUBLOSS(OBJ) returns resubstitution classification error for model
        %   OBJ.
        %
        %   L=RESUBLOSS(OBJ,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'LossFun'              - Function handle for loss, or string
        %                                'classiferror' representing a built-in
        %                                loss function. If you pass a function
        %                                handle FUN, RESUBLOSS calls it as shown
        %                                below:
        %                                      FUN(C,S,W,COST)
        %                                where C is an N-by-K logical matrix for N
        %                                rows in X and K classes in the ClassNames
        %                                property, S is an N-by-K numeric matrix, W
        %                                is a numeric vector with N elements, and
        %                                COST is a K-by-K numeric matrix. C has one
        %                                true per row for the true class. S is a
        %                                matrix of negated loss values for classes
        %                                with one row per observation, similar to
        %                                NEGLOSS output from PREDICT. W is a vector
        %                                of observation weights. COST is a matrix
        %                                of misclassification costs. Default:
        %                                'classiferror'
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
        %                                FUN, PREDICT calls it as shown below:
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
        %                                verbosity level, either 0 or 1. LOSS does
        %                                not display any diagnostic messages at
        %                                verbosity level 0 and displays diagnostic
        %                                messages at verbosity level 1. Default: 0
        %
        %   See also ClassificationECOC, loss.
            
            [varargout{1:nargout}] = ...
                resubLoss@classreg.learning.classif.FullClassificationModel(this,varargin{:});
        end
        
        function m = resubMargin(this,varargin)
        %RESUBMARGIN Classification margins by resubstitution.
        %   M=RESUBMARGIN(OBJ) returns resubstitution classification margins for
        %   model OBJ. Classification margin is the difference between
        %   classification score for the true class and maximal classification
        %   score for the false classes. The returned M is a numeric column-vector
        %   of length size(OBJ.X,1).
        %
        %   M=RESUBMARGIN(MODEL,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
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
        %                                FUN, PREDICT calls it as shown below:
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
        %                                RESUBMARGIN does not display any
        %                                diagnostic messages at verbosity level 0
        %                                and displays diagnostic messages at
        %                                verbosity level 1. Default: 0
        %
        %   See also ClassificationECOC, predict.
        
            m = resubMargin@classreg.learning.classif.FullClassificationModel(this,varargin{:});
        end
        
        function e = resubEdge(this,varargin)
        %RESUBEDGE Classification edge by resubstitution.
        %   E=RESUBEDGE(OBJ) returns classification edge obtained by model OBJ for
        %   training data OBJ.X and OBJ.Y. Classification edge is classification
        %   margin averaged over the entire data.
        %
        %   M=RESUBEDGE(MODEL,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
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
        %                                FUN, PREDICT calls it as shown below:
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
        %                                RESUBEDGE does not display any
        %                                diagnostic messages at verbosity level 0
        %                                and displays diagnostic messages at
        %                                verbosity level 1. Default: 0
        %
        %   See also ClassificationECOC, resubMargin, edge.
        
            e = resubEdge@classreg.learning.classif.FullClassificationModel(this,varargin{:});
        end        
    end
    
end


function [learners,weights] = localFitECOC(...
    tmp,M,X,C,W,prior,cost,doposterior,paropts,verbose)
%FITCECOC Fit an ECOC model by reducing to binary.
%   LEARNERS=LOCALFITECOC(TMP,M,X,C,W,PRIOR,COST,FITPOSTERIOR,PAROPTS,VERBOSE)
%   returns a cell array with L binary learners LEARNERS.
%       TMP     - Cell array of classification templates produced by a
%                 template function such as, for example, TEMPLATETREE with
%                 L elements. 
%       M       - K-by-L coding matrix for K classes and L binary
%                 classifiers 
%       X       - N-by-D predictor matrix 
%       C       - N-by-K logical matrix of class membership
%       W       - N-by-1 vector of observation weights
%       PRIOR   - 1-by-K vector of class prior probabilities
%       COST    - K-by-K matrix of misclassification costs
%       FITPOSTERIOR - Logical flag indicating if posterior probabilities
%                      need to be fitted.
%       PAROPTS - struct produced by STATSET to pass parallel computation
%                 options.
%       VERBOSE - Verbosity flag, 0 or 1.

[K,L] = size(M);

[N,K2] = size(C);
if K2~=K
    error(message('stats:ClassificationECOC:localFitECOC:MismatchBetweenCodingMatrixAndClassMembership',K,K2));
end

% Turn off SVM warning about perfect class separation for fitting posterior
% probabilities.
if doposterior
    f = @(z) strcmp(z.Method,'SVM');
    if any(cellfun(f,tmp))
        warning('off','stats:fitSVMPosterior:PerfectSeparation');
        cleanupObj = onCleanup(@() warning('on','stats:fitSVMPosterior:PerfectSeparation'));
    end
end

[useParallel,RNGscheme] = ...
    internal.stats.parallel.processParallelAndStreamOptions(paropts);

[learners,weights] = ...
    internal.stats.parallel.smartForSliceout(L,@loopBody,useParallel,RNGscheme);

weights = weights(:)';
    
    function [learner,weight] = loopBody(l,s)
        if isempty(s)
            s = RandStream.getGlobalStream;
        end
        
        pos = find(M(:,l)==+1);
        neg = find(M(:,l)==-1);
        
        if ~isempty(prior)
            lPrior = [sum(prior(neg)) sum(prior(pos))];
        else
            lPrior = 'empirical';
        end
        
        % If priors are 0 and 1, there is no point in applying the cost
        % correction. An edge case when a class has prior 1 and
        % misclassification cost 0 cannot be handled because only costs for
        % classes with non-zero priors are saved in the ClassSummary
        % property. When executing the crossval method, we cannot
        % reconstruct the non-zero costs passed to the model for
        % zero-probability classes.
        if ~isempty(cost) && all(lPrior>0)
            cost(isnan(cost)) = 0;
            lPrior = ...
                [sum(prior(neg))*prior(neg)*cost(neg,pos)*prior(pos)' ...
                 sum(prior(pos))*prior(pos)*cost(pos,neg)*prior(neg)'];            
        end
        
        idxpos = sum(C(:,pos),2)>0;
        idxneg = sum(C(:,neg),2)>0;
        
        y = zeros(N,1);
        y(idxpos) = +1;
        y(idxneg) = -1;
        
        x = X;
        w = W;
        x(y==0,:) = [];
        w(y==0)   = [];
        y(y==0)   = [];
        
        if verbose>0
            fprintf('%s\n',getString(message('stats:ClassificationECOC:localFitECOC:TrainingLearner',...
                l,tmp{l}.Method,L,sum(idxneg),sum(idxpos))));
            if verbose>1
                fprintf('%s',getString(message('stats:ClassificationECOC:localFitECOC:NegativeIndices')));
                for n=1:numel(neg)
                    fprintf(' %i',neg(n));
                end
                fprintf('\n');
                
                fprintf('%s',getString(message('stats:ClassificationECOC:localFitECOC:PositiveIndices')));
                for n=1:numel(pos)
                    fprintf(' %i',pos(n));
                end
                fprintf('\n\n');
            end
        end
        
        weight = sum(w);
        
        try
            full = fit(tmp{l},x,y,'ClassNames',[-1 1],'Prior',lPrior,'weights',w);
        catch me
            warning(message('stats:ClassificationECOC:localFitECOC:CannotFitLearner',...
                l,tmp{l}.Method,me.message));
            learner = {[]};
            return;
        end
        
        if doposterior && ~strcmp(full.ScoreType,'probability')
            if     isa(full,'ClassificationSVM')
                if verbose>0
                    fprintf('%s\n',getString(message('stats:ClassificationECOC:localFitECOC:FittingLearner',...
                        l,tmp{l}.Method)));
                end
                
                try
                    tb = tabulate(full.Y);
                    kfold = min(tb(:,2));
                    kfold = min(kfold,10);
                    if kfold<2
                        error(message('stats:ClassificationECOC:localFitECOC:NotEnoughDataToFitPosterior'));
                    end
                    cvpart = cvpartition(full.Y,'kfold',kfold,s);
                    full = fitPosterior(full,'cvpartition',cvpart);
                catch me
                    warning(message('stats:ClassificationECOC:localFitECOC:CannotFitPosterior',...
                        l,tmp{l}.Method,me.message));
                    learner = {[]};
                    return;
                end
                
            elseif isa(full,'ClassificationDiscriminant') ...
                    || isa(full,'ClassificationKNN') ...
                    || isa(full,'ClassificationTree')
                %             warning('Resetting ScoreTransform for learner %i (%s) for fitting posterior probabilities.',...
                %                 l,class(full));
                full.ScoreTransform = @classreg.learning.transform.identity;
                
            elseif isa(full,'classreg.learning.classif.ClassificationEnsemble') ...
                    && ~isempty(full.TransformToProbability)
                full.ScoreTransform = full.TransformToProbability;
                
                %                 && ismember(full.Method,{'AdaBoostM1' 'GentleBoost' 'LogitBoost' 'RUSBoost'})
                % %                 warning('Resetting ScoreTransform for learner %i (%s) for fitting posterior probabilities.',...
                % %                     l,class(full));
                %             full.ScoreTransform = @classreg.learning.transform.doublelogit;
                %
                %         elseif isa(full,'classreg.learning.classif.ClassificationEnsemble') ...
                %                 && ismember(full.Method,{'Bag' 'Subspace'})
                %             full.ScoreTransform = @classreg.learning.transform.identity;
                
            else
                error(message('stats:ClassificationECOC:localFitECOC:CannotUseLearnerToEstimatePosterior',...
                    l,tmp{l}.Method));
                
            end
        end
    
        learner = {compact(full)};
    end

end
