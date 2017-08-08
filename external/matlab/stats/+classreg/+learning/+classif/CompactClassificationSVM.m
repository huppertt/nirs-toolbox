classdef CompactClassificationSVM < classreg.learning.classif.ClassificationModel
%CompactClassificationSVM Support Vector Machine model for classification.
%   CompactClassificationSVM is an SVM model for classification with one or
%   two classes. This model can predict response for new data.
%
%   CompactClassificationSVM properties:
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
%
%   CompactClassificationSVM methods:
%       compareHoldout        - Compare two models using test data.
%       discardSupportVectors - Discard support vectors for linear SVM.
%       edge                  - Classification edge.
%       fitPosterior          - Find transformation from SVM scores to class posterior probabilities.
%       loss                  - Classification loss.
%       margin                - Classification margins.
%       predict               - Predicted response of this model.
%
%   See also ClassificationSVM.

%   Copyright 2013-2014 The MathWorks, Inc.

    properties(SetAccess=protected,GetAccess=public,Dependent=true)
        %ALPHA Coefficients obtained by solving the dual problem.
        %   The Alpha property is a vector with M positive elements for M support
        %   vectors. The PREDICT method computes scores F for predictor matrix X
        %   using
        %
        %   F = G(X,SupportVectors)*(SupportVectorLabels.*Alpha) + Bias
        %
        %   where G(X,SupportVectors) is an N-by-M matrix of kernel products for N
        %   rows in X and M rows in SupportVectors.
        %
        %   See also classreg.learning.classif.CompactClassificationSVM,
        %   Bias, SupportVectors, SupportVectorLabels, predict.
        Alpha;
        
        %BETA Coefficients for the primal linear problem.
        %   The Beta property is a vector of linear coefficients with P elements for P
        %   predictors. If this SVM model is obtained using a kernel function other than
        %   'linear', this property is empty. Otherwise scores F obtained by this model
        %   for predictor matrix X can be computed using
        %
        %   F = (X/S)*Beta + Bias
        %
        %   where S is the kernel scale saved in the KernelParameters.Scale property.
        %
        %   See also classreg.learning.classif.CompactClassificationSVM, Alpha, Bias,
        %   KernelParameters, SupportVectors, SupportVectorLabels, predict.
        Beta;
        
        %BIAS Bias term.
        %   The Bias property is a scalar specifying the bias term in the SVM
        %   model. The PREDICT method computes scores F for predictor matrix X
        %   using
        %
        %   F = G(X,SupportVectors)*(SupportVectorLabels.*Alpha) + Bias
        %
        %   where G(X,SupportVectors) is an N-by-M matrix of kernel products for N
        %   rows in X and M rows in SupportVectors.
        %
        %   See also classreg.learning.classif.CompactClassificationSVM, Alpha,
        %   Bias, SupportVectors, SupportVectorLabels, predict.
        Bias;
        
        %KERNELPARAMETERS Parameters of the kernel function.
        %   The KernelParameters is a struct with two fields:
        %       Function     - Name of the kernel function, a string.
        %       Scale        - Scale factor used to divide predictor values.
        %   The PREDICT method computes a kernel product between vectors x and z
        %   using Function(x/Scale,z/Scale).
        %
        %   See also classreg.learning.classif.CompactClassificationSVM, predict.
        KernelParameters;
        
        %MU Predictor means.
        %   The Mu property is either empty or a vector with P elements, one for
        %   each predictor. If training data were standardized, the Mu property is
        %   filled with means of predictors used for training. Otherwise the Mu
        %   property is empty. If the Mu property is not empty, the PREDICT method
        %   centers predictor matrix X by subtracting the respective element of Mu
        %   from every column.
        %
        %   See also classreg.learning.classif.CompactClassificationSVM, predict.
        Mu;
        
        %SIGMA Predictor standard deviations.
        %   The Sigma property is either empty or a vector with P elements, one for
        %   each predictor. If training data were standardized, the Sigma property
        %   is filled with standard deviations of predictors used for training.
        %   Otherwise the Sigma property is empty. If the Sigma property is not
        %   empty, the PREDICT method scales predictor matrix X by dividing every
        %   column by the respective element of Sigma (after centering).
        %
        %   See also classreg.learning.classif.CompactClassificationSVM, predict.
        Sigma;
        
        %SUPPORTVECTORS Support vectors.
        %   The SupportVectors property is an M-by-P matrix for M support vectors
        %   and P predictors. The PREDICT method computes scores F for predictor
        %   matrix X using
        %
        %   F = G(X,SupportVectors)*(SupportVectorLabels.*Alpha) + Bias
        %
        %   where G(X,SupportVectors) is an N-by-M matrix of kernel products for N
        %   rows in X and M rows in SupportVectors.
        %
        %   See also classreg.learning.classif.CompactClassificationSVM,
        %   Alpha, Bias, SupportVectorLabels, predict.
        SupportVectors;
        
        %SUPPORTVECTORLABELS Support vector labels.
        %   The SupportVectorLabels property is a vector with M elements for M
        %   support vectors saved in SupportVectors. For two-class learning with
        %   non-zero prior probabilities, an element of this vector is set to -1 if
        %   the respective observation is of class ClassNames(1) and +1 if the
        %   respective observation is of class ClassNames(2) in the training data.
        %   For one-class learning, all elements of this vector are set to +1. The
        %   PREDICT method computes scores F for predictor matrix X using
        %
        %   F = G(X,SupportVectors)*(SupportVectorLabels.*Alpha) + Bias
        %
        %   where G(X,SupportVectors) is an N-by-M matrix of kernel products for N
        %   rows in X and M rows in SupportVectors.
        %
        %   See also classreg.learning.classif.CompactClassificationSVM,
        %   ClassNames, Alpha, Bias, SupportVectors, predict.
        SupportVectorLabels;
    end
    
    methods
        function a = get.Alpha(this)
            a = this.Impl.Alpha;
        end
        
        function a = get.Bias(this)
            a = this.Impl.Bias;
        end
        
        function b = get.Beta(this)
            b = this.Impl.Beta;
        end
        
        function p = get.KernelParameters(this)
            p.Function = this.Impl.KernelParameters.Function;
            p.Scale    = this.Impl.KernelParameters.Scale;
            % Offset has been incorporated in the Bias already.
            % For the user, offset for prediction is zero.
            
            if     strcmpi(p.Function,'polynomial')
                p.Order = this.Impl.KernelParameters.PolyOrder;
            elseif strcmpi(p.Function,'sigmoid')
                p.Sigmoid = this.Impl.KernelParameters.Sigmoid;
            end
        end
        
        function a = get.Mu(this)
            a = this.Impl.Mu;
        end
        
        function a = get.Sigma(this)
            a = this.Impl.Sigma;
        end
        
        function a = get.SupportVectors(this)
            a = this.Impl.SupportVectors;
        end
        
        function a = get.SupportVectorLabels(this)
            a = this.Impl.SupportVectorLabels;
        end
    end
    
    methods
        function [varargout] = predict(this,X,varargin)
        %PREDICT Predict response of the SVM model.
        %   [LABEL,SCORE]=PREDICT(MODEL,X) returns predicted class labels and
        %   scores for SVM model MODEL and a matrix of predictors X. X must be a
        %   numeric matrix of size N-by-P, where P is the number of predictors used
        %   for training this model. Classification labels LABEL have the same type
        %   as Y used for training. Scores SCORE are an N-by-K numeric matrix for N
        %   observations and K classes. The predicted label is assigned to the
        %   class with the largest score.
        %
        %   See also classreg.learning.classif.CompactClassificationSVM.
            
            [varargout{1:nargout}] = ...
                predict@classreg.learning.classif.ClassificationModel(this,X,varargin{:});
        end
        
        function [varargout] = loss(this,X,Y,varargin)
        %LOSS Classification error.
        %   L=LOSS(MODEL,X,Y) returns classification error for SVM model MODEL
        %   computed using matrix of predictors X and true class labels Y. X must
        %   be a numeric matrix of size N-by-P, where P is the number of predictors
        %   used for training this model. Y must be of the same type as
        %   MODEL.ClassNames and have N elements.
        %
        %   L=LOSS(MODEL,X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'lossfun'          - Function handle for loss, or string
        %                            representing a built-in loss function.
        %                            Available loss functions for classification:
        %                            'binodeviance', 'classiferror', 'hinge', and
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
        %       'weights'         -  Vector of observation weights. By default the
        %                            weight of every observation is set to 1. The
        %                            length of this vector must be equal to the
        %                            number of rows in X.
        %
        %   See also classreg.learning.classif.CompactClassificationSVM,
        %   ClassNames, predict.

            [varargout{1:nargout}] = ...
                loss@classreg.learning.classif.ClassificationModel(this,X,Y,varargin{:});
        end
        
        function [obj,trans] = fitPosterior(obj,X,Y)
        %FITPOSTERIOR Fit posterior probabilities
        %   OBJ=FITPOSTERIOR(OBJ,X,Y) finds the optimal sigmoid transformation for
        %   SVM model OBJ using predictor matrix X and class labels Y. Pass X as a
        %   floating-point matrix with N rows. Pass Y as a character matrix with N
        %   rows or a vector with N elements of one of the following types: single,
        %   double, logical, categorical, and cell array of strings. FITPOSTERIOR
        %   returns an object of type CompactClassificationSVM and sets the
        %   ScoreTransform property of this object to the optimal transformation.
        %   The PREDICT method for the updated model then returns posterior
        %   probabilities for 2nd output.
        %
        %   [OBJ,TRANS]=FITPOSTERIOR(OBJ,X,Y) also returns TRANS, a struct with
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
        %   See also classreg.learning.classif.CompactClassificationSVM,
        %   ClassNames, ScoreTransform.

            [obj,trans] = fitSVMPosterior(obj,X,Y);
        end
        
        function this = discardSupportVectors(this)
        %DISCARDSUPPORTVECTORS Discard support vectors for linear SVM.
        %   OBJ=DISCARDSUPPORTVECTORS(OBJ) empties the Alpha, SupportVectors, and
        %   SupportVectorLabels properties of a linear SVM model. After these properties
        %   are emptied, the PREDICT method can compute SVM scores and labels using the
        %   Beta coefficients. DISCARDSUPPORTVECTORS errors for a non-linear SVM model.
        %
        %   See also Alpha, Beta, SupportVectors, SupportVectorLabels, predict.
        
            this.Impl = discardSupportVectors(this.Impl);
        end
    end
    
    methods(Access=protected)
        function cl = getContinuousLoss(this)
            cl = [];
            if     isequal(this.PrivScoreTransform,@classreg.learning.transform.identity)
                cl = @classreg.learning.loss.hinge;
            elseif strcmp(this.ScoreType,'probability')
                cl = @classreg.learning.loss.quadratic;
            end
        end
        
        function this = CompactClassificationSVM(...
                dataSummary,classSummary,scoreTransform,scoreType,impl)
            this = this@classreg.learning.classif.ClassificationModel(...
                dataSummary,classSummary,scoreTransform,scoreType);
            this.Impl = impl;
            this.DefaultLoss = @classreg.learning.loss.classiferror;
            this.LabelPredictor = @classreg.learning.classif.ClassificationModel.maxScore;
            this.DefaultScoreType = 'inf';
        end
        
        function S = score(this,X,varargin)
            f = score(this.Impl,X,varargin{:});
            
            classnames = this.ClassSummary.ClassNames;
            
            % Initialize all scores to minus the score for the positive
            % class and fill out the column for the positive class later.
            % When two classes are passed to fitcsvm, but only one of them
            % has non-zero probability, scores obtained by one-class
            % learning can be used for two-class predictions.
            S = repmat(-f,1,numel(classnames));

            % For one-class learning, the obtained score is for the only
            % class with non-zero probability. For two-class learning, the
            % obtained score is for the 2nd class (+1) with non-zero
            % probability.
            [~,loc] = ismember(this.ClassSummary.NonzeroProbClasses,classnames);
            if numel(loc)==1 % one-class learning
                S(:,loc) = f;
            else             % binary SVM
                S(:,loc(2)) = f;
            end
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.classif.ClassificationModel(this,s);
            hasAlpha = ~isempty(this.Alpha);
            
            if ~hasAlpha
                s.Beta                = this.Beta;
            end
            
            if hasAlpha
                s.Alpha               = this.Alpha;
            end
            
            s.Bias                    = this.Bias;
            s.KernelParameters        = this.KernelParameters;
            
            if ~isempty(this.Mu)
                s.Mu                  = this.Mu;
            end
            if ~isempty(this.Sigma)
                s.Sigma               = this.Sigma;
            end
            
            if hasAlpha
                s.SupportVectors      = this.SupportVectors;
            end
            
            if hasAlpha && numel(this.ClassSummary.NonzeroProbClasses)>1
                s.SupportVectorLabels = this.SupportVectorLabels;
            end
        end
    end
   
end
