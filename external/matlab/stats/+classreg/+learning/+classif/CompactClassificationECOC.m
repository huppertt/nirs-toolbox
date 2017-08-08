classdef CompactClassificationECOC < classreg.learning.classif.ClassificationModel
%CompactClassificationECOC Error-correcting output code (ECOC) model.
%   CompactClassificationECOC is an ECOC model for multiclass learning.
%   This model can predict response for new data.
%
%   CompactClassificationECOC properties:
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
%
%   CompactClassificationECOC methods:
%       compareHoldout        - Compare two models using test data.
%       discardSupportVectors - Discard support vectors for linear SVM.
%       edge                  - Classification edge.
%       loss                  - Classification loss.
%       margin                - Classification margins.
%       predict               - Predicted response of this model.
%
%   See also ClassificationECOC.

%   Copyright 2014 The MathWorks, Inc.
    
    properties(GetAccess=public,SetAccess=protected)
        %BINARYLEARNERS Binary learners.
        %   The BinaryLearners property is a cell array of trained binary learners.
        %   Element I of this array is the learner trained to solve the binary
        %   problem specified by column I of the coding matrix.
        %
        %   See also classreg.learning.classif.CompactClassificationECOC,
        %   CodingMatrix.
        BinaryLearners = {};
        
        %BINARYLOSS Default binary loss function for prediction.
        %   The BinaryLoss property is a string specifying the default function for
        %   computing loss incurred by each binary learner.
        %
        %   See also classreg.learning.classif.CompactClassificationECOC, predict.
        BinaryLoss = [];
        
        %CODINGMATRIX Coding matrix.
        %   The CodingMatrix property is a K-by-L matrix for K classes and L binary
        %   learners. Its elements take values -1, 0 or +1. If element (I,J) of
        %   this matrix is -1, class I is included in the negative class for binary
        %   learner J; if this element is +1, class I is included in the positive
        %   class for binary learner J; and if this element is 0, class I is not
        %   used for training binary learner J.
        %
        %   See also classreg.learning.classif.CompactClassificationECOC,
        %   BinaryLearners.
        CodingMatrix = [];
        
        %LEARNERWEIGHTS Weights for binary learners.
        %   The LearnerWeights property is a row-vector with L elements for L
        %   binary learners. Element I of this vector is the total weight of
        %   observations used to train binary learner I.
        %
        %   See also classreg.learning.classif.CompactClassificationECOC,
        %   BinaryLearners.
        LearnerWeights = {};
    end
    
    methods(Access=protected)
        function this = setScoreType(~,~) %#ok<STOUT>
            error(message('stats:classreg:learning:classif:CompactClassificationECOC:setScoreType:DoNotUseScoreType'));
        end
        
        function cl = getContinuousLoss(this) %#ok<STOUT,MANU>
            % ContinuousLoss is replaced by BinaryLoss
            error(message('stats:classreg:learning:classif:CompactClassificationECOC:getContinuousLoss:DoNotUseContinuousLoss'));
        end
        
        function this = CompactClassificationECOC(...
                dataSummary,classSummary,scoreTransform,learners,weights,M)
            this = this@classreg.learning.classif.ClassificationModel(...
                dataSummary,classSummary,scoreTransform,[]);
            this.LabelPredictor = @classreg.learning.classif.ClassificationModel.maxScore;
            this.DefaultLoss = @classreg.learning.loss.classiferror;
            this.BinaryLearners = learners;
            this.LearnerWeights = weights;
            this.CodingMatrix = M;
            [this.BinaryLoss,this.DefaultScoreType] = ...
                classreg.learning.classif.CompactClassificationECOC.analyzeLearners(learners);
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.classif.ClassificationModel(this,s);
            s.BinaryLearners = this.BinaryLearners;
            s.CodingMatrix = this.CodingMatrix;
        end
        
        function [labels,negloss,pscore,posterior] = predictEmptyX(this,X)
            D = numel(this.PredictorNames);
            if size(X,2)~=D
                error(message('stats:classreg:learning:classif:ClassificationModel:predictEmptyX:XSizeMismatch', D));
            end
            labels = repmat(this.ClassNames(1,:),0,1);
            K = numel(this.ClassSummary.ClassNames);
            negloss = NaN(0,K);
            pscore = NaN(0,K);
            posterior = NaN(0,K);
        end
        
        function [labels,negloss,pscore,posterior] = predictForEmptyLearners(this,X)
            D = numel(this.PredictorNames);
            if size(X,2)~=D
                error(message('stats:classreg:learning:classif:ClassificationModel:predictEmptyX:XSizeMismatch', D));
            end
            K = numel(this.ClassSummary.ClassNames);
            L = numel(this.BinaryLearners);
            N = size(X,1);
            
            [~,cls] = max(this.Prior);
            labels = repmat(this.ClassNames(cls,:),N,1);
            negloss = NaN(N,K);
            pscore = NaN(N,L);
            posterior = zeros(N,K);
            posterior(:,cls) = 1;
        end
        
        function [negloss,pscore] = score(...
                this,X,dist,ignorezeros,useParallel,verbose)
            % Init
            trained = this.BinaryLearners;
            M = this.CodingMatrix;
            
            % Ignore zeros in the coding matrix?
            if ignorezeros
                M(M==0) = NaN;
            end
            
            % Compute scores for positive class only
            pscore = localScore(X,trained,useParallel,verbose);
            
            if verbose>0
                fprintf('%s\n',getString(message('stats:classreg:learning:classif:CompactClassificationECOC:score:PredictionsComputed')));
            end
                        
            % Loss per observation
            negloss = -localLoss(dist,M,pscore,useParallel);
            if verbose>0
                fprintf('%s\n',getString(message('stats:classreg:learning:classif:CompactClassificationECOC:score:LossComputed')));
            end            
        end    
    end
    
        
    methods
        function [labels,negloss,pscore,posterior] = predict(this,X,varargin)
        %PREDICT Predict response of the ECOC model.
        %   [LABEL,NEGLOSS]=PREDICT(MODEL,X) returns predicted class labels LABEL
        %   and negated values of the average binary loss per class NEGLOSS for
        %   ECOC model MODEL and matrix of predictors X. X must be a numeric matrix
        %   of size N-by-P, where P is the number of predictors used for training
        %   this model. Classification labels LABEL have the same type as Y used
        %   for training. Negative loss values NEGLOSS are an N-by-K matrix for N
        %   observations and K classes. The predicted label is assigned to the
        %   class with the largest negated average binary loss, or equivalently
        %   smallest average binary loss.
        %
        %   [~,~,PBSCORE]=PREDICT(MODEL,X) also returns positive-class scores
        %   predicted by the binary learners, an N-by-L matrix for N observations
        %   and L binary learners.
        %
        %   [~,~,~,POSTERIOR]=PREDICT(MODEL,X) also returns posterior probability
        %   estimates, an N-by-K matrix for N observations and K classes. PREDICT
        %   cannot compute these estimates unless you passed 'FitPosterior' as true
        %   to FITCECOC. If you set 'FitPosterior' to false for FITCECOC and
        %   request 4th output, PREDICT throws an error.
        %
        %   [...]=PREDICT(MODEL,X,'PARAM1',val1,'PARAM2',val2,...) specifies
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
        %                                the Kullback-Leibler divergence. PREDICT
        %                                ignores this parameter unless you request
        %                                4th output and set 'PosteriorMethod' to
        %                                'kl'. Default: 0
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
        %                                learners. PREDICT ignores this parameter
        %                                unless you request 4th output. Default:
        %                                'kl'
        %       'Verbose'              - Non-negative integer specifying the
        %                                verbosity level, either 0 or 1. PREDICT
        %                                does not display any diagnostic messages
        %                                at verbosity level 0 and displays
        %                                diagnostic messages at verbosity level 1.
        %                                Default: 0
        %
        %   See also classreg.learning.classif.CompactClassificationECOC, fitcecoc,
        %   BinaryLoss, CodingMatrix, statset.            
            
            % Empty data
            if isempty(X)
                [labels,negloss,pscore,posterior] = predictEmptyX(this,X);
                return;
            end
                        
            % If all binary learners are empty, predict into the majority
            % class.
            if all(cellfun(@isempty,this.BinaryLearners))
                [labels,negloss,pscore,posterior] = predictForEmptyLearners(this,X);
                return;
            end

            % Get args
            args = {   'binaryloss'     'decoding' 'verbose' ...
                'posteriormethod' 'numklinitializations'           'options'};
            defs = {this.BinaryLoss 'lossweighted'         0 ...
                             'kl'                      0 statset('parallel')};
            [userloss,decoding,verbose,postmethod,numfits,paropts] = ...
                internal.stats.parseArgs(args,defs,varargin{:});

            % Can compute posterior probabilities?
            doposterior = nargout>3;
            if doposterior 
                if ~strcmp(this.ScoreType,'probability')
                    error(message('stats:classreg:learning:classif:CompactClassificationECOC:predict:CannotFitProbabilities'));
                end
            end
            
            % Use weighted averaging?
            if ~ischar(decoding)
                error(message('stats:classreg:learning:classif:CompactClassificationECOC:predict:BadDecodingType'));
            end
            allowedVals = {'LossBased' 'LossWeighted'};
            tf = strncmpi(decoding,allowedVals,length(decoding));
            if sum(tf)~=1
                error(message('stats:classreg:learning:classif:CompactClassificationECOC:predict:BadDecodingValue'));
            end
            ignorezeros = tf(2);
            
            % Use QP to fit probabilities?
            if doposterior
                if ~ischar(postmethod)
                    error(message('stats:classreg:learning:classif:CompactClassificationECOC:predict:BadPosteriorMethodType'));
                end
                allowedVals = {'QP' 'KL'};
                tf = strncmpi(postmethod,allowedVals,length(postmethod));
                if sum(tf)~=1
                    error(message('stats:classreg:learning:classif:CompactClassificationECOC:predict:BadPosteriorMethodValue'));
                end
                doquadprog = tf(1);
            end
            
            % How many extra attempts for KL fitting?
            if doposterior
                if ~isempty(numfits) && ...
                        (~isscalar(numfits) || ~isnumeric(numfits) ...
                        || numfits~=round(numfits) || numfits<0)
                    error(message('stats:classreg:learning:classif:CompactClassificationECOC:predict:BadNumKLInitializationsType'));
                end
                if doquadprog && numfits>0
                    error(message('stats:classreg:learning:classif:CompactClassificationECOC:predict:BadNumKLInitializationsValue'));
                end
            end
            
            % Make sure loss makes sense
            if ~isa(userloss,'function_handle')
                if isempty(this.BinaryLoss)
                    error(message('stats:classreg:learning:classif:CompactClassificationECOC:predict:MustProvideCustomBinaryLoss'));
                end
                
                allowedVals = {'hamming' 'linear' 'quadratic' 'exponential' 'binodeviance' 'hinge'};
                tf = strncmpi(userloss,allowedVals,length(userloss));
                if sum(tf)~=1
                    error(message('stats:classreg:learning:classif:CompactClassificationECOC:predict:BinaryLoss'));
                end
                userloss = allowedVals{tf};
                
                if strcmp(userloss,'quadratic') && ...
                        ~(strcmp(this.ScoreType,'01') || strcmp(this.ScoreType,'probability'))
                    error(message('stats:classreg:learning:classif:CompactClassificationECOC:predict:QuadraticLossForInfRange'));
                end
                
                if ismember(userloss,{'linear' 'exponential' 'binodeviance' 'hinge'}) ...
                        && ~strcmp(this.ScoreType,'inf')
                    error(message('stats:classreg:learning:classif:CompactClassificationECOC:predict:BadBinaryLossFor01Range',userloss));
                end
            end
            
            % Make a function for loss
            if isa(userloss,'function_handle')
                dist = userloss;
            else
                % M is a matrix and f is a row-vector
                switch userloss
                    case 'hamming'
                        switch this.ScoreType
                            case 'inf'
                                dist = @(M,f) nanmean( 1 - sign(bsxfun(@times,M,f)), 2 )/2;
                            case {'01' 'probability'}
                                dist = @(M,f) nanmean( 1 - sign(bsxfun(@times,M,2*f-1)), 2)/2;
                        end
                    case 'linear'       % range must be 'inf'
                        dist = @(M,f) nanmean( 1 - bsxfun(@times,M,f), 2 )/2;
                    case 'quadratic'    % range must be '01' or 'probability'
                        dist = @(M,f) nanmean( (1 - bsxfun(@times,M,2*f-1) ).^2, 2 )/2;
                    case 'exponential'  % range must be 'inf'
                        dist = @(M,f) nanmean( exp( -bsxfun(@times,M,f) ), 2 )/2;
                    case 'binodeviance' % range must be 'inf'
                        dist = @(M,f) nanmean(log( 1 + exp(-2*bsxfun(@times,M,f)) ),2)/(2*log(2));
                    case 'hinge'        % range must be 'inf'
                        dist = @(M,f) nanmean( max(0, 1-bsxfun(@times,M,f) ), 2 )/2;
                end
            end
            
            % Get parallel options
            [useParallel,RNGscheme] = ...
                internal.stats.parallel.processParallelAndStreamOptions(paropts);
            
            % Get loss and score for the positive class
            [negloss,pscore] = score(...
                this,X,dist,ignorezeros,useParallel,verbose);
            
            % Get class labels
            labels = this.LabelPredictor(this.ClassNames,...
                this.Prior,this.Cost,negloss,@(x)x);
            
            % Fit posterior probabilities
            if doposterior
                posterior = posteriorFromRatio(...
                    this.CodingMatrix,pscore,this.LearnerWeights,...
                    verbose,doquadprog,numfits,useParallel,RNGscheme);
            end
        end
        
        function m = margin(this,X,Y,varargin)
        %MARGIN Classification margins.
        %   M=MARGIN(MODEL,X,Y) returns classification margins obtained by MODEL
        %   for matrix of predictors X and class labels Y. X must be a numeric
        %   matrix of size N-by-P, where P is the number of predictors used for
        %   training this model. Y must be of the same type as MODEL.ClassNames and
        %   have N elements. Classification margin is the difference between
        %   classification score for the true class and maximal classification
        %   score for the false classes. The returned M is a numeric column-vector
        %   of length size(X,1).
        %
        %   M=MARGIN(MODEL,X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
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
        %                                verbosity level, either 0 or 1. MARGIN
        %                                does not display any diagnostic messages
        %                                at verbosity level 0 and displays
        %                                diagnostic messages at verbosity level 1.
        %                                Default: 0
        %
        %   See also classreg.learning.classif.CompactClassificationECOC,
        %   predict.
        
            m = margin@classreg.learning.classif.ClassificationModel(this,X,Y,varargin{:});
        end
        
        function e = edge(this,X,Y,varargin)
        %EDGE Classification edge.
        %   E=EDGE(MODEL,X,Y) returns classification edge obtained by MODEL for
        %   matrix of predictors X and class labels Y. X must be a numeric matrix
        %   of size N-by-P, where P is the number of predictors used for training
        %   this model. Y must be of the same type as MODEL.ClassNames and have N
        %   elements. Classification edge is classification margin averaged over
        %   the entire data.
        %
        %   E=EDGE(OBJ,X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'Weights'              - Observation weights, a numeric vector of length
        %                                size(X,1). By default, all observation
        %                                weights are set to 1. If you supply
        %                                weights, EDGE computes weighted
        %                                classification edge.
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
        %                                verbosity level, either 0 or 1. EDGE does
        %                                not display any diagnostic messages at
        %                                verbosity level 0 and displays diagnostic
        %                                messages at verbosity level 1. Default: 0
        %
        %   See also classreg.learning.classif.CompactClassificationECOC,
        %   classreg.learning.classif.CompactClassificationECOC/margin.
        
            e = edge@classreg.learning.classif.ClassificationModel(this,X,Y,varargin{:});
        end
        
        function l = loss(this,X,Y,varargin)
        %LOSS Classification error.
        %   L=LOSS(MODEL,X,Y) returns classification error for model MODEL computed
        %   using matrix of predictors X and true class labels Y. X must be a
        %   numeric matrix of size N-by-P, where P is the number of predictors used
        %   for training this model. Y must be of the same type as MODEL.ClassNames
        %   and have N elements.
        %
        %   L=LOSS(MODEL,X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'LossFun'              - Function handle for loss, or string
        %                                'classiferror' representing a built-in
        %                                loss function. If you pass a function
        %                                handle FUN, LOSS calls it as shown below:
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
        %       'Weights'              - Vector of observation weights. By default
        %                                the weight of every observation is set to
        %                                1. The length of this vector must be equal
        %                                to the number of rows in X.
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
        %   See also classreg.learning.classif.CompactClassificationECOC, 
        %   predict.
        
            l = loss@classreg.learning.classif.ClassificationModel(this,X,Y,varargin{:});
        end
        
        function this = discardSupportVectors(this)
        %DISCARDSUPPORTVECTORS Discard support vectors for linear SVM binary learners.
        %   OBJ=DISCARDSUPPORTVECTORS(OBJ) empties the Alpha, SupportVectors, and
        %   SupportVectorLabels properties of binary learners that are linear SVM
        %   models. After these properties are emptied, the PREDICT method can compute
        %   predictions from those binary learners using their Beta properties.
        %
        %   See also classreg.learning.classif.CompactClassificationSVM, BinaryLearners,
        %   predict.
        
            % Find linear SVM learners
            f = @(z) isa(z,'classreg.learning.classif.CompactClassificationSVM') ...
                && strcmp(z.KernelParameters.Function,'linear');
            isLinearSVM = cellfun(f,this.BinaryLearners);
            
            % If no such learners, warn
            if ~any(isLinearSVM)
                warning(message('stats:classreg:learning:classif:CompactClassificationECOC:discardSupportVectors:NoLinearSVMLearners'));
                return;
            end
        
            % Discard SV's for linear SVM's
            idxLinearSVM = find(isLinearSVM);
            for i=1:numel(idxLinearSVM)
                n = idxLinearSVM(i);
                this.BinaryLearners{n} = discardSupportVectors(this.BinaryLearners{n});
            end
        end
    end
    
    methods(Static=true,Hidden=true)
        function [lossType,scoreType] = analyzeLearners(learners)            
            if isempty(learners)
                lossType = 'hamming';
                scoreType = 'unknown';
                return;
            end

            % Fill score and object type per learner
            L = numel(learners);
            scoreTypes = repmat({''},L,1);
            lossTypes  = repmat({''},L,1);
            for l=1:L
                lrn = learners{l};
                if ~isempty(lrn)
                    scoreTypes(l) = {lrn.ScoreType};
                    lossTypes(l)  = {lossToString(lrn.ContinuousLoss)};
                end
            end
               
            % Figure out the score type
            if     ismember('unknown',scoreTypes)
                scoreType = 'unknown';
                warning(message('stats:classreg:learning:classif:CompactClassificationECOC:analyzeLearners:UnknownScoreType'));
            elseif (ismember('01',scoreTypes) || ismember('probability',scoreTypes)) ...
                    && ismember('inf',scoreTypes)
                scoreType = 'unknown';
                warning(message('stats:classreg:learning:classif:CompactClassificationECOC:analyzeLearners:ScoreRangeMismatch'));
            elseif ismember('01',scoreTypes)
                scoreType = '01';
            elseif ismember('probability',scoreTypes)
                scoreType = 'probability';
            elseif ismember('inf',scoreTypes)
                scoreType = 'inf';
            else
                % We get here if all learners are empty, that is, none was
                % trained successfully.
                scoreType = 'unknown';
            end

            %
            % Figure out the default loss.
            %            
            lossTypes(strcmp(lossTypes,'')) = [];
            lossTypes = unique(lossTypes);
            if     isempty(lossTypes) % if all learners are empty, no loss is set
                lossType = '';
                warning(message('stats:classreg:learning:classif:CompactClassificationECOC:analyzeLearners:AllLearnersEmpty'));
            elseif strcmp(scoreType,'unknown') % if scoreType is unknown, warning has been thrown already
                lossType = '';
            elseif numel(lossTypes)==1 % set the unique loss
                lossType = lossTypes{1};
            else % if there is no unique loss and scoreType is known, set to hamming
                lossType = 'hamming';
                warning(message('stats:classreg:learning:classif:CompactClassificationECOC:analyzeLearners:HammingLoss'));
            end
        end
    end
    
end


function str = lossToString(fhandle)
if     isequal(fhandle,@classreg.learning.loss.quadratic)
    str = 'quadratic';
elseif isequal(fhandle,@classreg.learning.loss.hinge)
    str = 'hinge';
elseif isequal(fhandle,@classreg.learning.loss.exponential)
    str = 'exponential';
elseif isequal(fhandle,@classreg.learning.loss.binodeviance)
    str = 'binodeviance';
else
    str = 'unknown';
end
end


function pscore = localScore(X,trained,useParallel,verbose)

N = size(X,1);
L = numel(trained);

%pscore = NaN(N,L);

pscore = ...
    internal.stats.parallel.smartForSliceout(L,@loopBody,useParallel);

%  For one observation, the positive score for a binary learner is scalar.
%  smartForSliceout by default concatenates scalars into a column-vector.
%  We want a row-vector.
if N==1
    pscore = pscore(:)';
end

    function lscore = loopBody(l,~)
        if verbose>1
            fprintf('%s\n',getString(message('stats:classreg:learning:classif:CompactClassificationECOC:localScore:ProcessingLearner',l)));
        end
        
        if isempty(trained{l})
            lscore = NaN(N,1);
        else
            [~,s] = predict(trained{l},X);
            lscore = s(:,2);
        end
    end
end


function vloss = localLoss(dist,M,pscore,useParallel)

N = size(pscore,1);

% dist ensures that vloss is a concatenation of row-vectors.
vloss = internal.stats.parallel.smartForSliceout(N,@loopBody,useParallel);

    function lloss = loopBody(n,~)
        lloss = dist(M,pscore(n,:))';
    end
end


function Phat = posteriorFromRatio(M,R,W,verbose,doquadprog,T,useParallel,RNGscheme)
%POSTERIORFROMRATIO Compute posterior probabilities given their ratios.
%   P=POSTERIORFROMRATIO(M,R) returns an N-by-K matrix of estimated class
%   posterior probabilities for indicator matrix M of size K-by-L for K
%   classes and L classifiers and matrix of posterior ratios R of size
%   N-by-L for N observations. Element R(n,l) gives the ratio of the
%   posterior probability of the positive class at observation n for
%   learner l over the total posterior probability of the positive and
%   negative classes at observation n for learner l. For observation n,
%   posterior probabilities P(n,:) are estimated by solving a QP problem.
%   You must pass M as a matrix filled with 0, +1 and -1.
%
%   P=POSTERIORFROMRATIO(M,R,W,VERBOSE,DOQP,USEPARALLEL) accepts
%       W           - Row-vector of total data weights for learners with L
%                     elements.
%       VERBOSE     - Verbosity flag, 0 or 1.
%       DOQP        - Logical flag. If true, solve the least-squares
%                     problem by QUADPROG; if false, minimize the
%                     Kullback-Leibler divergence.
%       USEPARALLEL - Logical flag. If true, use parallel computing.
%       RNGSCHEME   - RNG scheme

if nargin<3
    W = ones(1,size(R,2));
end

if nargin<4
    verbose = 0;
end

if nargin<5
    doquadprog = false;
end

if doquadprog && isempty(ver('Optim'))
    error(message('stats:classreg:learning:classif:CompactClassificationECOC:posteriorFromRatio:NeedOptim'));
end

K = size(M,1);
N = size(R,1);

Mminus        = M;
Mminus(M~=-1) = 0;
Mplus         = M;
Mplus(M~=+1)  = 0;

if verbose>0
    fprintf('%s\n',getString(message('stats:classreg:learning:classif:CompactClassificationECOC:posteriorFromRatio:ComputingPosteriorProbs')));
end

    function p = loopBodyQP(n,~)        
        p = NaN(1,K);
        
        r = R(n,:);
        
        igood = ~isnan(r);
        if ~any(igood)
            return;
        end
        
        Q = bsxfun(@times,Mminus(:,igood),r(igood)) + ...
            bsxfun(@times,Mplus(:,igood),1-r(igood));
        H = Q*Q'; % K-by-K
        [p,~,exitflag] = quadprog(H,zeros(K,1),[],[],ones(1,K),1,zeros(K,1),ones(K,1),[],opts);
        
        if exitflag~=1
            warning(message('stats:classreg:learning:classif:CompactClassificationECOC:posteriorFromRatio:QuadprogFails',n));
        end
        
        p = p';
    end


    function p = loopBodyKL(n,s)
        if isempty(s)
            s = RandStream.getGlobalStream;
        end

        p = NaN(1,K);
        
        r = R(n,:);
        
        igood = ~isnan(r);
        if ~any(igood)
            return;
        end

        phat = zeros(T+2,K);
        dist = zeros(T+2,1);
        p0 = rand(s,T,K);
        p0 = bsxfun(@rdivide,p0,sum(p0,2));
        
        % Random initialization
        for t=1:T
            [phat(t,:),dist(t)] = ...
                minimizeKL(r(igood),Mminus(:,igood),Mplus(:,igood),W(igood),p0(t,:)');
        end
        
        % Uniform initial estimates
        [phat(T+1,:),dist(T+1)] = ...
            minimizeKL(r(igood),Mminus(:,igood),Mplus(:,igood),W(igood),repmat(1/K,K,1));

        % Built-in estimates based on approximate lsqnonneg solution
        [phat(T+2,:),dist(T+2)] = ...
            minimizeKL(r(igood),Mminus(:,igood),Mplus(:,igood),W(igood));
        
        % Best solution
        [~,tmin] = min(dist);
        p = phat(tmin,:);
    end

if doquadprog % Solve the QP problem
    opts = optimoptions(@quadprog,...
        'Algorithm','interior-point-convex','Display','off');
    
    Phat = internal.stats.parallel.smartForSliceout(N,@loopBodyQP,useParallel);
        
else          % Minimize KL divergence
    
    Phat = internal.stats.parallel.smartForSliceout(N,@loopBodyKL,useParallel,RNGscheme);
    
end

end


function [p,dist] = minimizeKL(r,Mminus,Mplus,W,p0)

if nargin<5
    K = size(Mminus,1);
    
    M = Mminus + Mplus;
    M(M==-1) = 0;
    p = lsqnonneg(M',r');
    
    doquit = false;
    if     all(p==0)
        p = repmat(1/K,K,1);
        doquit = true;
    elseif sum(p>0)==1
        p(p>0) = 1;
        doquit = true;
    end
    
    if doquit
        rhat = sum(bsxfun(@times,Mplus,p));
        rhat = rhat ./ (rhat - sum(bsxfun(@times,Mminus,p)));        
        dist = KLdistance(r,rhat,W);
        return;
    end
    
    p = max(p,100*eps);
    p = p/sum(p);
    p(p>1) = 1;
else
    p = p0;
end

rhat = sum(bsxfun(@times,Mplus,p));
rhat = rhat ./ (rhat - sum(bsxfun(@times,Mminus,p)));

dist = KLdistance(r,rhat,W);

delta = Inf;

iter = 1;

while delta>1e-6 && iter<=1000
    iter = iter + 1;
    
    numer = sum(    bsxfun(@times,Mplus,W.*r) -    bsxfun(@times,Mminus,W.*(1-r)), 2 );
    denom = sum( bsxfun(@times,Mplus,W.*rhat) - bsxfun(@times,Mminus,W.*(1-rhat)), 2 );
    
    i = denom<=0 & numer>0;
    if any(i)
        p(i) = 1;
        p(~i) = 0;
    else
        j = denom>0;
        p(j) = p(j) .* numer(j)./denom(j);
        p(~j) = 0;
    end
    
    p = max(p,100*eps);
    p = p/sum(p);
    
    rhat = sum(bsxfun(@times,Mplus,p));
    rhat = rhat ./ (rhat - sum(bsxfun(@times,Mminus,p)));
    
    distnew = KLdistance(r,rhat,W);
    
    delta = dist - distnew;
    
    dist = distnew;
end

end


function dist = KLdistance(r,rhat,w)

i = r   > 100*eps;
dist = sum( w(i).* r(i).*log(r(i)./rhat(i)) );

i = 1-r > 100*eps;
dist = dist + sum( w(i).*(1-r(i)).*log((1-r(i))./(1-rhat(i))) );

end
