classdef CompactClassificationDiscriminant < classreg.learning.classif.ClassificationModel
%CompactClassificationDiscriminant Discriminant analysis.
%   CompactClassificationDiscriminant is a discriminant analysis model. It
%   can predict response for new data.
%
%   CompactClassificationDiscriminant properties:
%       PredictorNames        - Names of predictors used for this discriminant.
%       ResponseName          - Name of the response variable.
%       ClassNames            - Names of classes in Y.
%       Cost                  - Misclassification costs.
%       Prior                 - Prior class probabilities.
%       ScoreTransform        - Transformation applied to predicted classification scores.
%       DiscrimType           - Discriminant type.
%       Gamma                 - Regularization parameter for correlation matrix of predictors.
%       Delta                 - Threshold for linear coefficients.
%       MinGamma              - Minimal allowed Gamma.
%       DeltaPredictor        - Threshold at which predictor is eliminated from the model.
%       Mu                    - Matrix of class means.
%       Sigma                 - Within-class covariance matrix.
%       BetweenSigma          - Between-class covariance matrix.
%       LogDetSigma           - Natural log of determinant of within-class covariance matrix.
%       Coeffs                - Discriminant coefficients.
%
%   CompactClassificationDiscriminant methods:
%       compareHoldout        - Compare two models using test data.
%       edge                  - Classification edge.
%       loss                  - Classification loss.
%       logP                  - Log of the unconditional probability.
%       margin                - Classification margins.
%       mahal                 - Mahalanobis distance.
%       nLinearCoeffs         - Number of non-zero linear coefficients.
%       predict               - Predicted response of this discriminant.
%
%   See also ClassificationDiscriminant.

%   Copyright 2010-2014 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=public,Dependent=true)
        %DISCRIMTYPE Discriminant type.
        %   The DiscrimType property is a string with possible values 'linear',
        %   'pseudolinear', 'diaglinear', 'quadratic', 'pseudoquadratic' and
        %   'diagquadratic'. Assign to this property to change the type of
        %   discriminant analysis. You can switch from one linear type to another
        %   linear type, or you can switch from one quadratic type to another
        %   quadratic type. You cannot switch from a linear type to a quadratic
        %   type or from a quadratic type to a linear type.
        %
        %   See also classreg.learning.classif.CompactClassificationDiscriminant.
        DiscrimType;

        %GAMMA Regularization parameter for correlation matrix of predictors.
        %   The Gamma property is a scalar in the range between 0 and 1. The value
        %   of Gamma is added to the main diagonal of the correlation matrix of
        %   predictors to make this matrix non-singular. Assign to this property to
        %   change the amount of regularization. If you assign 1 for linear
        %   discriminant, the discriminant sets its type to 'diagLinear'. If you
        %   assign a value between MinGamma and 1 for linear discriminant, the
        %   discriminant sets its type to 'linear'. You cannot assign values below
        %   the value of the MinGamma property. For quadratic discriminant, you can
        %   assign either 0 (for DiscrimType 'quadratic') or 1 (for DiscrimType
        %   'diagQuadratic').
        %
        %   See also classreg.learning.classif.CompactClassificationDiscriminant,
        %   DiscrimType, Delta, MinGamma.
        Gamma;
        
        %DELTA Threshold for linear coefficients.
        %   The Delta property is a non-negative scalar. It defines a threshold on
        %   linear discriminant coefficients per class. If a magnitude of a linear
        %   coefficient is below this threshold, this coefficient is set to zero
        %   and the respective predictor can be dropped from the model. Assign a
        %   large value to this property to eliminate many predictors from the
        %   model. This property must be 0 for quadratic discriminant.
        %
        %   See also classreg.learning.classif.CompactClassificationDiscriminant,
        %   DiscrimType, Gamma, DeltaPredictor.
        Delta;
    end
    
    properties(GetAccess=public,SetAccess=protected)
        %COEFFS Discriminant coefficients.
        %   The Coeffs property is a K-by-K array of structs for K classes. Element
        %   (I,J) of this array holds coefficients for discrimination of class I
        %   against class J. Every element is a struct with fields
        %       DiscrimType      - Discriminant type, a string.
        %       Const            - The constant term.
        %       Linear           - A column-vector with P linear coefficients for P
        %                          predictors.
        %       Quadratic        - A P-by-P matrix with quadratic coefficients for
        %                          P predictors if DiscrimType is 'quadratic' or
        %                          'pseudoquadratic'. A 1-by-P vector if
        %                          DiscrimType is 'diagquadratic'. If DiscrimType
        %                          is not one of the quadratic types, this property
        %                          is absent.
        %       Class1           - Name of class I. This field is of the same type
        %                          as ClassNames.
        %       Class2           - Name of class J. This field is of the same type
        %                          as ClassNames.
        %
        %   See also classreg.learning.classif.CompactClassificationDiscriminant,
        %   DiscrimType, ClassNames.
        Coeffs = [];
    end
    
    properties(GetAccess=protected,SetAccess=protected)
        PrivBetweenSigma = [];
    end
    
    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        %MU Matrix of class means.
        %   The Mu property is a K-by-P matrix of class means for K classes and P
        %   predictors.
        %
        %   See also classreg.learning.classif.CompactClassificationDiscriminant.
        Mu = [];
        
        %SIGMA Within-class covariance matrix.
        %   The Sigma property is a P-by-P matrix for linear or pseudo-linear
        %   discriminant and a P-by-P-by-K array for quadratic or pseudo-quadratic
        %   discriminant for P predictors and K classes. If Sigma is P-by-P, it
        %   holds the pooled-in covariance matrix. If Sigma is P-by-P-by-K,
        %   Sigma(:,:,I) returns the covariance matrix for class I. For diagonal
        %   linear discriminant, Sigma is a row-vector with P elements holding
        %   variance for each predictor. For diagonal quadratic discriminant, Sigma
        %   is a 1-by-P-by-K array in which Sigma(1,I,J) is the variance of
        %   predictor I for class J.
        %
        %   To compute these matrices, FITCDISCR subtracts the class means from the
        %   input matrix of predictors and uses unbiased estimates of the
        %   covariance.
        %
        %   See also classreg.learning.classif.CompactClassificationDiscriminant,
        %   Mu, fitcdiscr.
        Sigma;
        
        %BETWEENSIGMA Between-class covariance matrix.
        %   The BetweenSigma property is a P-by-P matrix with an unbiased estimate
        %   of the between-class covariance. BetweenSigma is computed from class
        %   means stored in the Mu property and class weights passed to FITCDISCR.
        %
        %   See also classreg.learning.classif.CompactClassificationDiscriminant,
        %   Mu, fitcdiscr.
        BetweenSigma;
        
        %LOGDETSIGMA Natural log of determinant of within-class covariance matrix.
        %   The LogDetSigma property is a numeric scalar for linear discriminant
        %   and a column-vector with K elements for quadratic discriminant
        %   with K classes.
        %
        %   See also classreg.learning.classif.CompactClassificationDiscriminant,
        %   Sigma.
        LogDetSigma;
        
        %MINGAMMA Minimal allowed value of Gamma.
        %   The MinGamma property is a scalar between 0 and 1. It is the minimal
        %   value of Gamma required for inverting the correlation matrix of
        %   predictors.
        %
        %   See also classreg.learning.classif.CompactClassificationDiscriminant,
        %   Gamma.
        MinGamma;
        
        %DELTAPREDICTOR Threshold at which predictor is eliminated from the
        %model.
        %   The DeltaPredictor property is a row-vector with P elements for P
        %   predictors. If an element of this vector is below the value of the
        %   Delta property, the linear coefficient for this predictor is zero. For
        %   quadratic discriminant, this property is a vector with P zeros.
        %
        %   See also classreg.learning.classif.CompactClassificationDiscriminant,
        %   Delta.
        DeltaPredictor;
    end
    
    methods(Static=true,Access=protected)
        function [s,v,d] = decompose(sigma)
            d = diag(sigma)';
            p = numel(d);
            if any( d<0 & abs(d)>p*eps(max(d)) )
                error(message('stats:classreg:learning:classif:CompactClassificationDiscriminant:decompose:NegativeD'));
            end
            d(d<0) = 0;
            d = sqrt(d);
            badD = d<p*eps(max(d));
            d(badD) = 0;
            invD = 1./d;
            invD(badD) = 0;
            sigma = bsxfun(@times,invD',bsxfun(@times,sigma,invD));
            sigma = (sigma+sigma')/2;
            [v,s] = eig(sigma);
            s = diag(s);
            if any( s<0 & abs(s)>p*eps(max(s)) )
                error(message('stats:classreg:learning:classif:CompactClassificationDiscriminant:decompose:NegativeS'));
            end
            s(s<0) = 0;
            s = sqrt(s);
            s(s<p*eps(max(s))) = 0;
        end
    end
    
    methods
        function mu = get.Mu(this)
            mu = this.Impl.Mu;
        end
        
        function dt = get.DiscrimType(this)
            dt = this.Impl.Type;
        end
        
        function this = set.DiscrimType(this,discrimType)
            this = setDiscrimType(this,discrimType);
        end
        
        function gamma = get.Gamma(this)
            gamma = this.Impl.Gamma;
        end
        
        function delta = get.Delta(this)
            delta = this.Impl.Delta;
        end
        
        function this = set.Gamma(this,gamma)
            this = setGamma(this,gamma);
        end
        
        function this = set.Delta(this,delta)
            this = setDelta(this,delta);
        end
        
        function logsig = get.LogDetSigma(this)
            logsig = this.Impl.LogDetSigma;
        end
                
        function sig = get.Sigma(this)
            sig = this.Impl.Sigma;
        end
        
        function S = get.BetweenSigma(this)
            % If BetweenSigma was passed by the user, return it
            if ~isempty(this.PrivBetweenSigma)
                S = this.PrivBetweenSigma;
                return;
            end
            
            % Otherwise compute BetweenSigma from class means. Means for
            % classes without observations have been set to NaN's during
            % training.
            Wj = this.Impl.ClassWeights;
            centeredMu = this.Impl.CenteredMu;
            weightedMu = bsxfun(@times,centeredMu,sqrt(Wj));
            Wj = Wj(~isnan(Wj));
            S = weightedMu'*weightedMu / (1-Wj'*Wj);            
        end
        
        function mingam = get.MinGamma(this)
            mingam = this.Impl.MinGamma;
        end
        
        function delpred = get.DeltaPredictor(this)
            delpred = deltaPredictor(this.Impl);
        end
    end
    
    methods(Hidden)
        % Provide these as public hidden methods to allow direct assignment
        % to cross-validated folds in cvshrink method of the full class        
        function this = setGamma(this,gamma)
            this = setBaseGamma(this,gamma);
        end
        
        function this = setDelta(this,delta,checkValidity)
            if nargin<3
                checkValidity = true;
            end
            this = setBaseDelta(this,delta,checkValidity);
        end        
    end
    
    methods
        function [label,scores,cost] = predict(this,X)
        %PREDICT Predict response of the discriminant.
        %   [LABEL,POSTERIOR,COST]=PREDICT(DISCR,X) returns predicted class labels
        %   LABEL, posterior probabilities POSTERIOR and misclassification costs
        %   COST for discriminant DISCR and a matrix of predictors X. X must be a
        %   floating-point matrix of size N-by-P, where P is the number of
        %   predictors used for training this model. Classification labels LABEL
        %   have the same type as Y used for training. Posterior probabilities
        %   POSTERIOR are an N-by-K numeric matrix for N observations and K
        %   classes. COST is an N-by-K matrix with predicted misclassification
        %   costs per class. The predicted label is assigned to the class with the
        %   minimal misclassification cost.
        %
        %   See also classreg.learning.classif.CompactClassificationDiscriminant,
        %   mahal.
  
            % Empty data
            if isempty(X)
                [label,scores,cost] = predictEmptyX(this,X);
                return;
            end

            % Get log(P(x,k)) from score() method
            logP = score(this,X);
            
            % Divide by the largest probability to avoid overflow
            maxLogP = max(logP,[],2);
            P = exp(bsxfun(@minus,logP,maxLogP));
            
            % Get posterior probabilities
            sumP = nansum(P,2); % P(x)
            posterior = bsxfun(@times,P,1./(sumP)); % P(k|x)
            
            % Transform posterior, compute expected cost and find the class
            % with minimal predicted cost
            [label,scores,cost] = this.LabelPredictor(this.ClassNames,...
                this.Prior,this.Cost,posterior,this.PrivScoreTransform);
         end
        
        function M = mahal(this,X,varargin)
        %MAHAL Mahalanobis distance to class means.
        %   M=MAHAL(DISCR,X) returns squared Mahalanobis distance from observations
        %   in X to the class means for discriminant DISCR. X must be a
        %   floating-point matrix of size N-by-P, where P is the number of
        %   predictors used for training this discriminant. M is a numeric matrix
        %   of size N-by-K for K classes.
        %
        %   M=MAHAL(DISCR,X,'classlabels',Y) returns a column-vector M with N
        %   elements. Element I in this vector is the squared Mahalanobis distance
        %   from the I-th row of X to the mean for class defined by the I-th
        %   element of Y. Y must be of the same type as DISCR.Y and have N
        %   elements.
        %
        %   See also classreg.learning.classif.CompactClassificationDiscriminant.
            
            if ~isfloat(X) || ~ismatrix(X)
                error(message('stats:classreg:learning:classif:CompactClassificationDiscriminant:mahal:BadX'));
            end
            
            [N,p] = size(X);
            
            if p~=numel(this.PredictorNames)
                error(message('stats:classreg:learning:classif:CompactClassificationDiscriminant:mahal:XPredictorMismatch', numel( this.PredictorNames )));
            end
            
            % Decode input args
            args = {'classlabels'};
            defs = {           []};
            Y = internal.stats.parseArgs(args,defs,varargin{:});
            
            % Convert labels into a class count matrix
            C = [];
            if ~isempty(Y)
                Y = classreg.learning.internal.ClassLabel(Y);
                C = classreg.learning.internal.classCount(this.ClassSummary.ClassNames,Y);
                if size(C,1)~=N
                    error(message('stats:classreg:learning:classif:CompactClassificationDiscriminant:mahal:MismatchSizeXandY', N));
                end
            end

            % Returns scores for all known classes but compute them only
            % for classes with observations in the training data
            K = length(this.ClassSummary.ClassNames);                                    
            [~,Knonempty] = ismember(this.ClassSummary.NonzeroProbClasses,...
                this.ClassSummary.ClassNames);
            Knonempty = Knonempty(:)';

            % For observations with NaN's, log-scores are -Inf and
            % predicted posteriors are 0.
            nanX = any(isnan(X),2);
            X(nanX,:) = [];
            mah = NaN(size(X,1),K);
                                
            % Loop over classes
            mah(:,Knonempty) = mahal(this.Impl,Knonempty,X);
            
            % Return computed Mahalanobis distance per class
            M = NaN(N,K);
            M(~nanX,:) = mah;
            if ~isempty(C)
                M = sum(M.*C,2);
            end
        end
        
        function logp = logP(this,X)
        %LOGP Log of the unconditional probability density.
        %   LOGP=LOGP(DISCR,X) returns a column-vector with N elements for N rows
        %   in matrix of predictors X. Element I in this vector is the log of the
        %   unconditional probability of observing row I of X computed using the
        %   discriminant model DISCR.
        %
        %   See also classreg.learning.classif.CompactClassificationDiscriminant,
        %   mahal.
            
            % Get log(P(x,k)) from score() method
            logP = score(this,X);
            
            % Divide by the largest probability to avoid overflow
            maxLogP = max(logP,[],2);
            P = exp(bsxfun(@minus,logP,maxLogP));
            
            % Take log of the sum of P(x,k) over k
            p = size(X,2);
            logp = log(nansum(P,2)) + maxLogP - .5*p*log(2*pi);
        end
        
        function nCoeffs = nLinearCoeffs(this,delta)
        %NLINEARCOEFFS Number of non-zero linear coefficients.
        %   NCOEFFS=NLINEARCOEFFS(DISCR) returns the number of non-zero linear
        %   coefficients in linear discriminant DISCR.
        %
        %   NCOEFFS=NLINEARCOEFFS(DISCR,DELTA) returns the number of non-zero
        %   linear coefficients in discriminant DISCR for threshold parameter
        %   DELTA.
        %
        %   For quadratic discriminant, NCOEFFS is always set to the total number
        %   of predictors.
        %
        %   See also classreg.learning.classif.CompactClassificationDiscriminant,
        %   Delta.
            
            if nargin<2
                delta = this.Delta;
            end
            nCoeffs = nLinearCoeffs(this.Impl,delta);
        end
    end
         
    methods(Access=protected)
        function this = CompactClassificationDiscriminant(...
                dataSummary,classSummary,scoreTransform,scoreType,trained)
            this = this@classreg.learning.classif.ClassificationModel(...
                dataSummary,classSummary,scoreTransform,scoreType);
            if ~isempty(trained)
                this.PrivBetweenSigma   = trained.BetweenSigma;
                this.Impl               = trained.Impl;
                if trained.FillCoeffs
                    this                = fillCoeffs(this);
                    this                = adjustConstTerms(this);
                end
            end
        end
            
        function logP = score(this,X,varargin)
            % Get Mahalanobis distances to all classes
            mah = mahal(this,X);
            
            % Get log(P(x,k)) from Mahalanobis distance
            prior = this.Prior;
            logp = bsxfun(@minus,log(prior)-this.LogDetSigma'/2,mah/2);
            logp(:,prior==0) = -Inf;
            
            % Return computed log(P(x,k)) per class. For observations with
            % NaN's, log-scores are -Inf and predicted posteriors are 0.
            % Mahalanobis distance for observations with NaN's are NaN's,
            % and they need to be mapped to -Inf for log(P(x,k)).
            [N,K] = size(mah);
            nanX = any(isnan(X),2);
            logP = -Inf(N,K);
            logP(~nanX,:) = logp(~nanX,:);
        end
        
        function this = setPrior(this,prior)
            this = setPrivatePrior(this,prior);
            this = adjustConstTerms(this);
        end
        
        function this = setCost(this,cost)
            this = setPrivateCost(this,cost);
            this = adjustConstTerms(this);            
        end
        
        function this = setBaseDiscrimType(this,discrimType)
            this.Impl = setType(this.Impl,discrimType);
        end

        % setDiscrimType method in the derived full class does the same
        % plus sets DiscrimType in ModelParams to allow cross-validation.
        function this = setDiscrimType(this,discrimType)
            this = setBaseDiscrimType(this,discrimType);
            if ~isempty(this.Coeffs)
                this = fillCoeffs(this);
                this = adjustConstTerms(this);
            end
        end
                
        function this = setBaseGamma(this,gamma)
            if ~isnumeric(gamma) || ~isscalar(gamma) || gamma<0 || gamma>1 || isnan(gamma)
               error(message('stats:classreg:learning:classif:CompactClassificationDiscriminant:setBaseGamma:BadGamma'));
            end
            this.Impl = setGamma(this.Impl,gamma);
            if ~isempty(this.Coeffs)
                this = fillCoeffs(this);
                this = adjustConstTerms(this);
            end
        end
        
        function this = setBaseDelta(this,delta,checkValidity)
            if ~isnumeric(delta) || ~isscalar(delta) || delta<0 || isnan(delta)
               error(message('stats:classreg:learning:classif:CompactClassificationDiscriminant:setBaseDelta:BadDelta'));
            end
            this.Impl = setDelta(this.Impl,delta);
            
            % checkValidity should be true when the user works with a
            % discriminant object and false for cross-validation, ensemble
            % learning etc - in all situations when repetitive warning
            % messages would be annoying.
            if checkValidity && delta>0
                delrange = deltaRange(this.Impl);
                if delta<delrange(1) || delta>delrange(2)
                    warning(message('stats:classreg:learning:classif:CompactClassificationDiscriminant:setBaseDelta:DeltaOutOfRange', sprintf( '%g', delrange( 1 ) ), sprintf( '%g', delrange( 2 ) )));
                end
            end
            if ~isempty(this.Coeffs)
                this = fillCoeffs(this);
                this = adjustConstTerms(this);
            end
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.classif.ClassificationModel(this,s);
            s.DiscrimType = this.DiscrimType;
            s.Mu          = this.Mu;
            s.Coeffs      = this.Coeffs;
        end

        % Adjust intercepts every time the user assigns into Prior or Cost
        % properties. The intercepts are adjusted only for non-empty
        % classes (classes with observations found in Y during training).
        %
        % Increasing Prior(k) or Cost(k,:) for class k leads to pushing the
        % separating hyperplane away from the mean for class k. This is
        % equivalent to increasing Coeffs(k,:).Const.
        function this = adjustConstTerms(this)
            % No need to adjust const terms if coeffs were not filled in
            % fitting.
            if isempty(this.Coeffs)
                return;
            end
            
            % Get all class pairs
            [~,Knonempty] = ismember(this.ClassSummary.NonzeroProbClasses,...
                this.ClassSummary.ClassNames);
            Knonempty = Knonempty(:)';
            pairs = combnk(Knonempty,2)';
            npairs = size(pairs,2);
            
            % Copy the new prior and cost
            prior = this.Prior;
            cost = this.Cost;
            
            % Loop over the class pairs and reset the const terms
            for k=1:npairs
                i = pairs(1,k);
                j = pairs(2,k);
                this.Coeffs(i,j).Const = -constantTerm(this.Impl,i,j) ...
                    + log(prior(i)) - log(prior(j)) ...
                    + log(cost(i,j)) - log(cost(j,i));
                this.Coeffs(j,i).Const = -this.Coeffs(i,j).Const;
            end
        end

        % Compute coefficients from class means and covariance matrices.
        % Computation is done for non-empty classes only (classes with
        % observations found in Y during training).
        function this = fillCoeffs(this)
            K = length(this.ClassSummary.ClassNames);
            p = numel(this.PredictorNames);
            [~,Knonempty] = ismember(this.ClassSummary.NonzeroProbClasses,...
                this.ClassSummary.ClassNames);
            Knonempty = Knonempty(:)';

            % LDA or QDA?
            isquadratic = ~isempty(strfind(lower(this.DiscrimType),'quadratic'));

            % Set up the default empty struct
            s = struct('DiscrimType',this.DiscrimType,'Const',[],'Linear',[]);            
            if isquadratic
                s.Quadratic = [];
            end
            coeffs = repmat(s,K,K);
            s.DiscrimType = '';
            coeffs(1:K+1:end) = s;

            % Get all class pairs
            pairs = combnk(Knonempty,2)';
            npairs = size(pairs,2);
            
            % Initialize terms
            C = zeros(1,npairs); % constant term
            L = zeros(p,npairs); % linear term
            if ~strcmpi(this.DiscrimType,'diagquadratic')
                Q = zeros(p,p,npairs);
            else
                Q = zeros(1,p,npairs);
            end
            
            for j=1:npairs
                i1 = pairs(1,j);
                i2 = pairs(2,j);
                C(j) = constantTerm(this.Impl,i1,i2);
                L(:,j) = linearCoeffs(this.Impl,i1,i2);
                if isquadratic
                    Q(:,:,j) = quadraticCoeffs(this.Impl,i1,i2);
                end
            end            

            % Copy the computed coeffs into the struct
            for k=1:npairs
                i = pairs(1,k);
                j = pairs(2,k);

                coeffs(i,j).Const = -C(k);
                coeffs(i,j).Linear = L(:,k);
                
                coeffs(j,i).Const = C(k);
                coeffs(j,i).Linear = -L(:,k);
                
                if isquadratic
                    coeffs(i,j).Quadratic = Q(:,:,k);
                    coeffs(j,i).Quadratic = -Q(:,:,k);
                end
            end

            % Copy class names into the struct. Do this even for empty
            % classes fr which coeffs are not computed.
            pairs = combnk(1:K,2)';
            npairs = size(pairs,2);
            classnames = this.ClassSummary.ClassNames;
            for k=1:npairs
                i = pairs(1,k);
                j = pairs(2,k);
                coeffs(i,j).Class1 = labels(classnames(i));
                coeffs(i,j).Class2 = labels(classnames(j));
                coeffs(j,i).Class1 = labels(classnames(j));
                coeffs(j,i).Class2 = labels(classnames(i));
            end            

            for k=1:K
                coeffs(k,k).Class1 = labels(classnames(k));
                coeffs(k,k).Class2 = labels(classnames(k));
            end

            % Store coeffs in the object
            this.Coeffs = coeffs;
        end
    end
    
    methods(Static,Hidden)
        function this = loadobj(obj)
            if isempty(obj.Impl)
                % Load 11b discriminant
                trained.BetweenSigma = obj.PrivBetweenSigma;
                trained.FillCoeffs = ~isempty(obj.Coeffs);
                isquadratic = ~isempty(strfind(lower(obj.PrivDiscrimType),...
                    'quadratic'));
                if isquadratic
                    trained.Impl = ...
                        classreg.learning.impl.QuadraticDiscriminantImpl(...
                        obj.PrivDiscrimType,obj.D,obj.S,obj.V,0,0,obj.Mu,...
                        ones(size(obj.Mu,1),1),false);
                else
                    trained.Impl = ...
                        classreg.learning.impl.LinearDiscriminantImpl(...
                        obj.PrivDiscrimType,obj.D,obj.S,obj.V,0,0,obj.Mu,...
                        ones(size(obj.Mu,1),1),false);
                end
                this = classreg.learning.classif.CompactClassificationDiscriminant(...
                    obj.DataSummary,obj.ClassSummary,obj.PrivScoreTransform,[],...
                    trained);
            else
                % Load post-11b discriminant
                this = obj;
            end
        end
    end
    
end
