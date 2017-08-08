classdef ClassificationDiscriminant < ...
        classreg.learning.classif.FullClassificationModel & ...
        classreg.learning.classif.CompactClassificationDiscriminant
%ClassificationDiscriminant Discriminant analysis.
%   ClassificationDiscriminant is a discriminant model for classification.
%   It can predict response for new data. It also stores data used for
%   training and can compute resubstitution predictions.
%
%   An object of this class cannot be created by calling the constructor.
%   Use FITCDISCR to create a ClassificationDiscriminant object by fitting
%   the discriminant to training data.
%
%   This class is derived from CompactClassificationDiscriminant.
%
%   ClassificationDiscriminant properties:
%       NumObservations       - Number of observations.
%       X                     - Matrix of predictors used to train this discriminant.
%       XCentered             - Matrix of predictors with subtracted class means.
%       Y                     - True class labels used to train this discriminant.
%       W                     - Weights of observations used to train this discriminant.
%       ModelParameters       - Discriminant parameters.
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
%   ClassificationDiscriminant methods:
%       compact               - Compact this discriminant.
%       compareHoldout        - Compare two models using test data.
%       crossval              - Cross-validate this discriminant.
%       cvshrink              - Cross-validate regularization of this discriminant.
%       edge                  - Classification edge.
%       logP                  - Log of the unconditional probability.
%       loss                  - Classification loss.
%       margin                - Classification margins.
%       mahal                 - Mahalanobis distance.
%       nLinearCoeffs         - Number of non-zero linear coefficients.
%       predict               - Predicted response of this discriminant.
%       resubEdge             - Resubstitution classification edge.
%       resubLoss             - Resubstitution classification loss.
%       resubMargin           - Resubstitution classification margins.
%       resubPredict          - Resubstitution predicted response.
%
%   Example: Grow a classification discriminant for Fisher's iris data.
%       load fisheriris
%       d = fitcdiscr(meas,species,'PredictorNames',{'SL' 'SW' 'PL' 'PW'})
%
%   See also fitcdiscr,
%   classreg.learning.classif.CompactClassificationDiscriminant.
    
%   Copyright 2011-2014 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        %XCENTERED Matrix of predictors with subtracted class means.
        %   The XCentered property is a numeric matrix of size N-by-P, where N is
        %   the number of observations (rows) and P is the number of predictors
        %   (columns) in the training data. I-th row of this matrix is obtained by
        %   subtracting the mean vector for class Y(I) from I-th row of X.
        %
        %   See also ClassificationDiscriminant, X, Y, Mu.
        XCentered;
    end
    
    methods
        function x = get.XCentered(this)
            gidx = grp2idx(this.PrivY,this.ClassSummary.ClassNames);
            x = this.X - this.Mu(gidx,:);
        end
    end
    
    methods(Static,Hidden)
        function temp = template(varargin)
            classreg.learning.FitTemplate.catchType(varargin{:});
            temp = classreg.learning.FitTemplate.make('Discriminant','type','classification',varargin{:});
        end
        
        function this = fit(X,Y,varargin)
            % Exclude prior and cost from fit arguments
            args = {'prior' 'cost'};
            defs = {     []     []};
            [prior,cost,~,fitArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});

            % Fit
            temp = classreg.learning.FitTemplate.make(...
                'Discriminant','type','classification',fitArgs{:});
            this = fit(temp,X,Y);
            
            % Set the correct label predictor for cross-validated model
            this.LabelPredictor = @classreg.learning.classif.ClassificationModel.minCost;
 
            % Process prior and cost
            if ~isempty(prior) && ~strncmpi(prior,'empirical',length(prior))
                this.Prior = prior;
            end
            if ~isempty(cost)
                this.Cost = cost;
            end
        end

        function this = make(mu,sigma,varargin)            
            % Check type of mu and sigma
            if ~isnumeric(mu) || ~isnumeric(sigma)
                error(message('stats:ClassificationDiscriminant:make:MuSigmaNotNumeric'));
            end
            
            % Check dimensionality of mu and sigma
            if ~ismatrix(mu) || ndims(sigma)>3
                error(message('stats:ClassificationDiscriminant:make:MuSigmaBadSize'));
            end
            
            % Verify dimensions of mu and sigma match
            [K,p] = size(mu);
            if ismatrix(sigma) 
                if ~all(size(sigma)==[p p])
                    error(message('stats:ClassificationDiscriminant:make:SizeLDASigmaMismatch'));
                end
            else
                if ~all(size(sigma)==[p p K])
                    error(message('stats:ClassificationDiscriminant:make:SizeQDASigmaMismatch'));
                end
            end
            
            % Check for NaN's
            if any(isnan(mu(:))) || any(isnan(sigma(:)))
                error(message('stats:ClassificationDiscriminant:make:MuSigmaWithNaNs'));
            end
            
            % Verify that the covariance matrix is symmetric
            if ismatrix(sigma)
                if ~all(all( abs(sigma-sigma') < p*eps(max(abs(diag(sigma)))) ))
                    error(message('stats:ClassificationDiscriminant:make:PooledInSigmaNotSymmetric'));
                end
            else
                for k=1:K
                    if ~all(all( abs(sigma(:,:,k)-sigma(:,:,k)') < p*eps(max(abs(diag(sigma(:,:,k))))) ))
                        error(message('stats:ClassificationDiscriminant:make:ClassSigmaNotSymmetric', k));
                    end
                end
            end
            
            % Check optional args
            args = {'betweensigma' 'classnames' 'predictornames' ...
                'responsename' 'prior' 'cost' 'fillcoeffs'};
            defs = {            []           []               {} ...
                            ''      []     []         'on'};
            [betweenSigma,classnames,predictornames,responsename,prior,cost,fillCoeffs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Verify between-class covariance matrix
            if ~isempty(betweenSigma)
                if ~isnumeric(betweenSigma) || ~ismatrix(betweenSigma) ...
                        || ~all(size(betweenSigma)==[p p])
                    error(message('stats:ClassificationDiscriminant:make:BadBetweenSigma', p, p));
                end
                if ~all(all( abs(betweenSigma-betweenSigma') < max(K,p)*eps(max(abs(diag(betweenSigma)))) ))
                    error(message('stats:ClassificationDiscriminant:make:BetweenSigmaNotSymmetric'));
                end
                betweenSigma = (betweenSigma+betweenSigma')/2;
                d = eig(betweenSigma);
                if any( d<0 & abs(d)>max(K,p)*eps(max(d)) )
                    error(message('stats:ClassificationDiscriminant:make:BetweenSigmaNotPosSemiDefinite'));
                end
                d( d<10*max(K,p)*eps(max(d)) ) = 0;
                r = sum(d>0); % rank of betweenSigma
                if r>K-1
                    error(message('stats:ClassificationDiscriminant:make:RankBetweenSigmaTooLarge', r, K - 1, K));
                end
            end
            
            % Process class names
            if isempty(classnames)
                classnames = 1:K;
            end
            classnames = classreg.learning.internal.ClassLabel(classnames);
            if numel(classnames)~=K
                error(message('stats:ClassificationDiscriminant:make:ClassNamesSizeMismatch', K));
            end
            
            % Process predictor names
            if isempty(predictornames)
                predictornames = classreg.learning.internal.defaultPredictorNames(p);
            else
                if ~iscellstr(predictornames)
                    if ~ischar(predictornames)
                        error(message('stats:ClassificationDiscriminant:make:BadPredictorType'));
                    end
                    predictornames = cellstr(predictornames);
                end
                if length(predictornames)~=p
                    error(message('stats:ClassificationDiscriminant:make:PredictorMismatch', p));
                end
            end
            predictornames = predictornames(:)';
            
            % Process response name
            if isempty(responsename)
                responsename = 'Y';
            else
                if ~ischar(responsename)
                    error(message('stats:ClassificationDiscriminant:make:BadResponseName'));
                end
            end
            
            % Process prior and cost
            if ~isempty(prior) && ischar(prior) ...
                    && strncmpi(prior,'empirical',length(prior))
               error(message('stats:ClassificationDiscriminant:make:EmpiricalPriorNotAllowed'));
            end
            prior = classreg.learning.classif.FullClassificationModel.processPrior(...
                prior,ones(1,K),classnames,classnames);
            prior = prior/sum(prior);
            cost = classreg.learning.classif.FullClassificationModel.processCost(...
                cost,prior,classnames,classnames);
            
            % Process FillCoeffs
            if ~ischar(fillCoeffs) || ...
                    ~(strcmpi(fillCoeffs,'on') || strcmpi(fillCoeffs,'off'))
                error(message('stats:ClassificationDiscriminant:make:BadFillCoeffs'));
            end
            fillCoeffs = strcmpi(fillCoeffs,'on');
            
            % Decompose the covariance matrix
            isquadratic = ~ismatrix(sigma);
            if isquadratic
                d = zeros(1,p,K);
                s = cell(K,1);
                v = cell(K,1);
                for k=1:K
                    [s{k},v{k},d(1,:,k)] = ...
                        classreg.learning.classif.CompactClassificationDiscriminant.decompose(sigma(:,:,k));
                end
                discrimType = 'quadratic';
                if any(d(:)==0)
                    discrimType = 'pseudoQuadratic';
                else
                    for k=1:K
                        if any(s{k}==0)
                            discrimType = 'pseudoQuadratic';
                        end
                    end
                end
            else
                [s,v,d] = classreg.learning.classif.CompactClassificationDiscriminant.decompose(sigma);
                discrimType = 'linear';
                if any(s==0) || any(d==0)
                    discrimType = 'pseudoLinear';
                end
            end
            
            % Make compact object
            dataSummary.PredictorNames          = predictornames;
            dataSummary.CategoricalPredictors   = [];
            dataSummary.ResponseName            = responsename;
            classSummary.ClassNames             = classnames;
            classSummary.NonzeroProbClasses     = classnames;
            classSummary.Prior                  = prior;
            classSummary.Cost                   = cost;
            scoreTransform                      = @classreg.learning.transform.identity;
            if isquadratic
                trained.Impl                    = ...
                    classreg.learning.impl.QuadraticDiscriminantImpl(...
                    discrimType,d,s,v,0,0,mu,ones(K,1),false);
            else
                trained.Impl                    = ...
                    classreg.learning.impl.LinearDiscriminantImpl(...
                    discrimType,d,s,v,0,0,mu,ones(K,1),false);
            end
            trained.BetweenSigma                = betweenSigma;
            trained.FillCoeffs                  = fillCoeffs;
            this = classreg.learning.classif.CompactClassificationDiscriminant(...
                dataSummary,classSummary,scoreTransform,[],trained);
        end
    end
    
    methods(Hidden)
        function this = ClassificationDiscriminant(X,Y,W,modelParams,...
                dataSummary,classSummary,scoreTransform)
            % Protect against being called by the user
            if nargin~=7 || ischar(W)
                error(message('stats:ClassificationDiscriminant:ClassificationDiscriminant:DoNotUseConstructor'));
            end
            
            % The checks on categorical predictors and floating-point X
            % should be imposed in prepareData method. At present
            % ClassificationDiscriminant has prepareData method inherited
            % from the base class. I do not want to introduce a new method
            % just to address these two checks. Will revise if the list
            % gets longer. inarsky 10/24/2012 
            
            % Categorical predictors not allowed
            if ~isempty(dataSummary.CategoricalPredictors)
                error(message('stats:ClassificationDiscriminant:ClassificationDiscriminant:CategPred'));
            end
            
            % Accept only floating-point data
            if ~isfloat(X)
                error(message('stats:ClassificationDiscriminant:ClassificationDiscriminant:BadXType'));
            end

            % Base constructor
            this = this@classreg.learning.classif.FullClassificationModel(...
                X,Y,W,modelParams,dataSummary,classSummary,scoreTransform);
            this = this@classreg.learning.classif.CompactClassificationDiscriminant(...
                dataSummary,classSummary,scoreTransform,[],[]);
            
            % Discriminant mode
            discrimType = this.ModelParams.DiscrimType;
            if isempty(strfind(lower(discrimType),'linear'))
                isquadratic = true;
            else
                isquadratic = false;
            end
            
            % For quadratic DA, all classes must have weights
            classHasWeight = ...
                ismember(this.ClassSummary.ClassNames,this.ClassSummary.NonzeroProbClasses);
            if isquadratic && ~all(classHasWeight)
                classname = cellstr(this.ClassSummary.ClassNames(find(~classHasWeight,1)));
                error(message('stats:ClassificationDiscriminant:ClassificationDiscriminant:QuadDAwithZeroWeightClassNotAllowed',classname{1}));
            end

            % Get class means
            gidx = grp2idx(this.PrivY,this.ClassSummary.NonzeroProbClasses);
            K = numel(this.ClassSummary.NonzeroProbClasses); % must be equal to max(gidx)
            p = size(this.X,2);
            gmeans = NaN(K,p);
            for k=1:K
                gmeans(k,:) = classreg.learning.internal.wnanmean(this.X(gidx==k,:),this.W(gidx==k));
                idxnan = find(isnan(gmeans(k,:)),1);
                if ~isempty(idxnan)
                    classname = cellstr(this.ClassSummary.NonzeroProbClasses(k));
                    error(message('stats:ClassificationDiscriminant:ClassificationDiscriminant:NaNPredictor', idxnan, classname{ 1 }));
                end
            end
                        
            % Get max likelihood estimates of the covariance matrix
            % diagonals
            W = this.W;
            if isquadratic
                D = zeros(1,p,K);
                invD = zeros(1,p,K);
                for k=1:K
                    d = sqrt(classreg.learning.internal.wnanvar(this.X(gidx==k,:),W(gidx==k),0));
                    badD = (d <= sum(gidx==k)*eps(max(d))) | isnan(d);
                    d(badD) = 0;
                    if strcmpi(discrimType,'quadratic') && any(d==0)
                        classname = cellstr(this.ClassSummary.ClassNames(k));
                        error(message('stats:ClassificationDiscriminant:ClassificationDiscriminant:ZeroDiagCovQuad', this.PredictorNames{ find( d==0, 1 ) }, classname{ 1 }));
                    end
                    D(1,:,k) = d;
                    invD(1,:,k) = 1./d;
                    invD(1,badD,k) = 0;
                end
            else
                D = sqrt(classreg.learning.internal.wnanvar(this.X-gmeans(gidx,:),W,0));
                badD = (D <= numel(W)*eps(max(D))) | isnan(D);
                D(badD) = 0;
                if strcmpi(discrimType,'linear') && any(D==0)
                    error(message('stats:ClassificationDiscriminant:ClassificationDiscriminant:ZeroDiagCovLin', this.PredictorNames{ find( D==0, 1 ) }));
                end                
                invD = 1./D;
                invD(badD) = 0;
            end
            
            % Remove all NaN's
            nanX = any(isnan(this.X),2);
            if any(nanX)
                this.X(nanX,:)   = [];
                this.PrivY(nanX) = [];
                this.W(nanX)     = [];
            end
            if isempty(this.X)
                error(message('stats:ClassificationDiscriminant:ClassificationDiscriminant:NoDataAfterNaNsRemoved'));
            end
                        
            % Renormalize weights
            this.W = this.W/sum(this.W);
            
            % Save initial values.
            nonzeroClassNames = this.ClassSummary.NonzeroProbClasses;
            prior = this.ClassSummary.Prior;
            cost = this.ClassSummary.Cost;
            if isempty(cost)
                K = numel(nonzeroClassNames);
                cost = ones(K) - eye(K);
            end
            
            % Get matrix of class weights
            C = classreg.learning.internal.classCount(nonzeroClassNames,this.PrivY);
            WC = bsxfun(@times,C,this.W);
            Wj = sum(WC,1);
            gmeans(Wj==0,:) = [];
            
            % Remove observations for classes with zero prior probabilities
            % or costs
            [this.X,this.PrivY,C,WC,Wj,prior,cost,nonzeroClassNames] = ...
                classreg.learning.classif.FullClassificationModel.removeZeroPriorAndCost(...
                this.X,this.PrivY,C,WC,Wj,prior,cost,nonzeroClassNames);
            
            % Normalize priors in such a way that the priors in present
            % classes add up to one.  Normalize weights to add up to the
            % prior in the respective class.
            prior = prior/sum(prior);
            WC = bsxfun(@times,WC,prior./Wj);
            this.W = sum(WC,2);
            
            % Compute weights and group indices
            W = this.W;
            Wj2 = sum(WC.^2,1);
            gidx = grp2idx(this.PrivY,nonzeroClassNames);
            K = numel(Wj);
            
            % Put processed values into summary structure
            this.ClassSummary = ...
                classreg.learning.classif.FullClassificationModel.makeClassSummary(...
                this.ClassSummary.ClassNames,nonzeroClassNames,prior,cost);
            
            % For pooled-in covariance estimate, check overall data size
            N = size(this.X,1);
            if N<=K && ~isquadratic
                error(message('stats:ClassificationDiscriminant:ClassificationDiscriminant:TooFewObsForLinear'));
            end
            
            % For quadratic discriminants, get the group size and
            % normalize observation weights within each group.
            if isquadratic
                gsize = sum(C,1);
                idx = find(gsize<2,1);
                if ~isempty(idx)
                    classname = cellstr(this.ClassSummary.NonzeroProbClasses(idx));
                    error(message('stats:ClassificationDiscriminant:ClassificationDiscriminant:TooFewObsForQuad', classname{ 1 }));
                end
                % The bias correction for the estimate of the covariance
                % matrix for class k is given by a sum of squared weights
                % in this class under constraint that sum of weights is 1.
                WC = bsxfun(@rdivide,WC,Wj);
                W = sum(WC,2);
                Wj2 = sum(WC.^2,1);
            end
            C = []; %#ok<NASGU>
            WC = []; %#ok<NASGU>
            
            % Get gamma and delta
            gamma = this.ModelParams.Gamma;
            delta = this.ModelParams.Delta;

            % Training amounts to centering X by class and then decomposing
            % into a product of 3 matrices: 
            % D=std(X); [~,S,V]=svd(X*diag(1./D)) 
            %
            % If the input gamma is 0 (default) and the matrix is singular,
            % Gamma property in the constructed object is set to MinGamma
            % but this value is not copied to ModelParams. Assume that in
            % this case the user wants to cross-validate or perform
            % ensemble learning by putting the minimal amount of
            % regularization into every discriminant. If the user passes
            % positive gamma below MinGamma, he gets an error.
            if isquadratic
                S = cell(K,1);
                V = cell(K,1);
                for k=1:K
                    [~,s,v] = svd(bsxfun(@times,...
                        bsxfun(@times,...
                        bsxfun(@minus,this.X(gidx==k,:),gmeans(k,:)),...
                        sqrt(W(gidx==k))),...
                        invD(1,:,k)), 'econ');
                    s = diag(s);
                    bads = abs(s)<max(gsize(k),p)*eps(max(abs(s)));
                    s(bads) = 0;
                    S{k} = s;
                    biasCorr = sqrt(1-Wj2(k));
                    D(1,:,k) = D(1,:,k)/biasCorr;
                    V{k} = v;
                end
            else
                [~,S,V] = svd(bsxfun(@times, bsxfun(@times,this.X-gmeans(gidx,:),sqrt(W)), invD), 'econ');
                S = diag(S);
                badS = abs(S)<max(N,p)*eps(max(abs(S)));
                S(badS) = 0;
                biasCorr = sqrt(1-sum(Wj2./Wj));
                D = D/biasCorr;
            end
                        
            % Compute class means and class weights.
            %
            % For DA, unlike for most classifiers, we do not rid of classes
            % with zero priors and misclassification costs.
            % NonzeroProbClasses list has a special meaning: It includes
            % classes with zero priors but excludes classes without
            % observations in training Y. The rationale is this:            
            %   - For LDA use classes with zero priors to compute the
            %        pooled-in covariance matrix.
            %   - For QDA compute covariance matrices for classes with zero
            %        priors and all the coefficients anyway.
            % Nothing can be ever predicted into a zero-prior class and the
            % const term separating a zero-prior class from another class
            % is Inf (or -Inf).
            %
            % NonzeroProbClasses is a subset of ClassNames (all classes
            % passed in by the user which may include classes without
            % observations in Y). For classes without observations the
            % class means are set to NaN's.
            mu = NaN(numel(this.ClassSummary.ClassNames),p);
            [~,pos] = ismember(this.ClassSummary.NonzeroProbClasses,this.ClassSummary.ClassNames);
            mu(pos,:) = gmeans;
            classWeights = zeros(numel(this.ClassSummary.ClassNames),1);
            classWeights(pos) = Wj';
            
            % Make the implementation object
            if isquadratic
                this.Impl = classreg.learning.impl.QuadraticDiscriminantImpl(...
                    discrimType,D,S,V,gamma,delta,mu,classWeights,this.ModelParams.SaveMemory);
            else
                this.Impl = classreg.learning.impl.LinearDiscriminantImpl(...
                    discrimType,D,S,V,gamma,delta,mu,classWeights,this.ModelParams.SaveMemory);
            end
            
            % Fill coefficients if desired
            if this.ModelParams.FillCoeffs
                this               = fillCoeffs(this);
                this               = adjustConstTerms(this);
            end
        end
        
        function this = setGamma(this,gamma)
            this = setBaseGamma(this,gamma);
            this.ModelParams.Gamma = this.Impl.Gamma;
            this.ModelParams.DiscrimType = this.Impl.Type;
        end
        
        function this = setDelta(this,delta,checkValidity)
            if nargin<3
                checkValidity = true;
            end
            this = setBaseDelta(this,delta,checkValidity);
            this.ModelParams.Delta = this.Impl.Delta;
            this.ModelParams.DiscrimType = this.Impl.Type;
        end        
    end
    
    methods(Access=protected)
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.classif.FullClassificationModel(this,s);
            s = propsForDisp@classreg.learning.classif.CompactClassificationDiscriminant(this,s);
        end
        
        function this = setDiscrimType(this,mode)
            % Almost like the parent method in the compact class, but also
            % set DiscrimType in ModelParams, so it can be used for
            % cross-validation.
            this = setBaseDiscrimType(this,mode);
            if ~isempty(this.Coeffs)
                this = fillCoeffs(this);
                this = adjustConstTerms(this);
            end
            this.ModelParams.Gamma = this.Impl.Gamma;
            this.ModelParams.Delta = this.Impl.Delta;
            this.ModelParams.DiscrimType = this.Impl.Type;
        end                
    end
    
    methods
        function partModel = crossval(this,varargin)
        %CROSSVAL Cross-validate this discriminant.
        %   CVDISCR=CROSSVAL(DISCR) builds a partitioned discriminant from
        %   discriminant DISCR represented by a full object for classification. You
        %   can then assess the predictive performance of this discriminant on
        %   cross-validated data using methods and properties of CVDISCR. By
        %   default, CVDISCR is built using 10-fold cross-validation on the
        %   training data. CVDISCR is of class ClassificationPartitionedModel.
        %
        %   CVDISCR=CROSSVAL(DISCR,'PARAM1',val1,'PARAM2',val2,...) specifies
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
        %      'SaveMemory'  - If 'on', defers computing the full covariance matrix
        %                      in every fold until it is needed for prediction. If
        %                      'off', computes and stores the full covariance
        %                      matrix for every fold. Set this parameter to 'on' if
        %                      X has thousands of predictors. Default: setting used
        %                      to obtain DISCR
        %
        %   See also ClassificationDiscriminant, cvpartition,
        %   classreg.learning.partition.ClassificationPartitionedModel.

            % Catch base arguments
            forbiddenArgs = [classreg.learning.FitTemplate.AllowedBaseFitObjectArgs ...
                {'discrimtype' 'fillcoeffs'}];
            idxBaseArg = find(ismember(lower(varargin(1:2:end)),forbiddenArgs));
            if ~isempty(idxBaseArg)
                error(message('stats:ClassificationDiscriminant:crossval:NoBaseArgs',...
                    varargin{2*idxBaseArg-1}));
            end
            
            % Save memory?
            args = {'savememory'};
            defs = {          ''};
            [savememory,~,extraArgs] = internal.stats.parseArgs(args,defs,varargin{:});
            if ~isempty(savememory)
                if ~ischar(savememory) || (~strcmpi(savememory,'on') && ~strcmpi(savememory,'off'))
                    error(message('stats:ClassificationDiscriminant:crossval:BadSaveMemory'));
                end
                savememory = strcmpi(savememory,'on');
            end
            
            % Fix ModelParams
            modelparams = this.ModelParams;
            modelparams.FillCoeffs = false;
            if ~isempty(savememory)
                modelparams.SaveMemory = savememory;
            end
            temp = classreg.learning.FitTemplate.make(this.ModelParams.Method,...
                'type','classification','scoretransform',this.PrivScoreTransform,...
                'modelparams',modelparams,'CrossVal','on',extraArgs{:});
            
            % Fit arguments do not include prior and cost because DA
            % training does not depend on these. They are set in the
            % trained object afterwards.
            partModel = fit(temp,this.X,this.Y,'Weights',this.W,...
                'predictornames',this.DataSummary.PredictorNames,...
                'responsename',this.ResponseName,'classnames',this.ClassNames);
            partModel.ScoreType = this.ScoreType;
            
            % Prior and cost cannot be passed to fit() because classes with
            % zero prior and zero cost (both zero column and zero row) are
            % removed by the ensemble. Pass the default prior and cost to
            % fit() and then assign afterwards.
            partModel.Prior = this.Prior;
            partModel.Cost = this.Cost;
        end
    
        function cmp = compact(this)
        %COMPACT Compact discriminant analysis.
        %   CMP=COMPACT(DISCR) returns an object of class
        %   CompactClassificationDiscriminant holding the structure of the trained
        %   discriminant. The compact object does not contain X and Y used for
        %   training.
        %
        %   See also ClassificationDiscriminant,
        %   classreg.learning.classif.CompactClassificationDiscriminant.

            trained = struct('Impl',this.Impl,...
                'BetweenSigma',[],'FillCoeffs',this.ModelParams.FillCoeffs);
            cmp = classreg.learning.classif.CompactClassificationDiscriminant(...
                this.DataSummary,this.ClassSummary,...
                this.PrivScoreTransform,this.PrivScoreType,trained);
        end

        function [mcr,gamma,delta,npred] = cvshrink(this,varargin)
        %CVSHRINK Cross-validate regularization of linear discriminant.
        %   ERR=CVSHRINK(DISCR) returns a vector of cross-validated classification
        %   error values for 11 values of the regularization parameter Gamma
        %   uniformly spaced between MinGamma and 1.
        %
        %   [ERR,GAMMA]=CVSHRINK(DISCR) also returns a vector of used Gamma values.
        %
        %   [ERR,GAMMA,DELTA,NUMPRED]=CVSHRINK(DISCR) also returns matrix DELTA and
        %   matrix NUMPRED. By default DELTA is a zero scalar and NUMPRED is a
        %   row-vector with 11 elements set to the number of predictors (columns)
        %   in property X.
        %
        %   If you pass optional arguments 'gamma' or 'NumGamma', ERR and NUMPRED
        %   are row-vectors of length L for L elements in GAMMA. If you pass
        %   optional arguments 'delta' or 'NumDelta', ERR, DELTA and NUMPRED are
        %   matrices of size L-by-M for L elements in GAMMA and M thresholds on
        %   linear coefficients. ERR(I,J) is cross-validated error for discriminant
        %   with Gamma set to GAMMA(I) and Delta set to DELTA(I,J). NUMPRED(I,J) is
        %   the number of predictors with non-zero linear coefficients for
        %   discriminant with Gamma set to GAMMA(I) and Delta set to DELTA(I,J).
        %
        %   CVSHRINK performs 10-fold cross-validation by default. Use one of the
        %   following options for different kinds of cross-validation:
        %   'CVPartition', 'Holdout', 'KFold', and 'Leaveout'.
        %
        %   ERR=CVSHRINK(DISCR,'PARAM1',val1,'PARAM2',val2,...) specifies
        %   optional parameter name/value pairs:
        %       'Gamma'       - Vector of Gamma values for cross-validation.
        %                       Default: 0:0.1:1
        %       'NumGamma'    - Number of Gamma intervals for cross-validation.
        %                       CVSHRINK cross-validates discriminant using
        %                       NumGamma+1 values of Gamma uniformly spaced from
        %                       MinGamma to 1. If you pass 'gamma', CVSHRINK ignores
        %                       'NumGamma'. Default: 10
        %       'Delta'       - Row-vector or matrix of Delta values for
        %                       cross-validation. If you pass 'delta' as a
        %                       row-vector with M elements, CVSHRINK returns matrix
        %                       DELTA with M columns and fills DELTA(:,J) with
        %                       element J of the input vector. If you pass 'delta'
        %                       as a matrix, the number of rows in this matrix must
        %                       be equal to the number of elements in 'gamma'. In
        %                       this case, CVSHRINK copies the input matrix into the
        %                       returned matrix DELTA. Default: 0
        %       'NumDelta'    - Number of Delta intervals for cross-validation.
        %                       For every value of Gamma, CVSHRINK cross-validates
        %                       discriminant using NumDelta+1 values of Delta
        %                       uniformly spaced from zero to the maximal Delta at
        %                       which all predictors are eliminated for this value
        %                       of Gamma. If you pass 'delta', CVSHRINK ignores
        %                       'NumDelta'. Default: 0
        %       'CVPartition' - A partition created with CVPARTITION to use in
        %                       cross-validated discriminant.
        %       'Holdout'     - Holdout validation uses the specified fraction of
        %                       the data for test, and uses the rest of the data
        %                       for training. Specify a numeric scalar between 0
        %                       and 1.
        %       'KFold'       - Number of folds to use in cross-validated
        %                       discriminant, a positive integer. Default: 10
        %       'Leaveout'    - Use leave-one-out cross-validation by setting to
        %                       'on'. Not recommended for data with many
        %                       predictors.
        %       'Verbose'     - 0, 1 or 2. Increase the verbosity level to see more
        %                       progress messages from CVSHRINK. Default: 0
        %
        %   Example:
        %       load ovariancancer
        %       d = fitcdiscr(obs,grp,'FillCoeffs','off','SaveMemory','on')
        %       [err,gamma,delta,npred] = cvshrink(d,'NumGamma',30,'NumDelta',30,'Verbose',1)
        %
        %       % Plot classification error
        %       figure
        %       imagesc(1:31,gamma,err)
        %       colorbar
        %       xlabel('Delta steps')
        %       ylabel('Gamma')
        %       title('Classification error')
        %
        %       % Plot number of predictors in the model
        %       figure
        %       imagesc(1:31,gamma,npred)
        %       colorbar
        %       xlabel('Delta steps')
        %       ylabel('Gamma')
        %       title('Number of predictors')
        %
        %   See also ClassificationDiscriminant,
        %   classreg.learning.classif.CompactClassificationDiscriminant, fit,
        %   Gamma, MinGamma, Delta, DeltaPredictor.
            
            % Get input args
            args = {'gamma' 'NumGamma' 'delta' 'NumDelta' 'Verbose'};
            defs = {     []         10       0         []         0};
            [gamma,ngamma,delta,ndelta,verbose,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});

            % Print out?
            if ~isnumeric(verbose) || ~isscalar(verbose) || verbose<0
                error(message('stats:ClassificationDiscriminant:cvshrink:BadPrint'));
            end
            verbose = ceil(verbose);
            
            % Get cross-validation options
            [~,partitionArgs,extraArgs] = ...
                classreg.learning.generator.Partitioner.processArgs(extraArgs{:},'CrossVal','on');
            if ~isempty(extraArgs)
                error(message('stats:ClassificationDiscriminant:cvshrink:TooManyOptionalArgs'));
            end
            
            C = classreg.learning.internal.classCount(this.ClassSummary.ClassNames,this.PrivY);
            cost = this.Cost;
            
            % Cross-validate
            cv = crossval(this,'savememory','on',partitionArgs{:});
            if verbose>0
                fprintf(1,'%s\n',getString(message('stats:ClassificationDiscriminant:cvshrink_DoneCrossvalidating')));
            end
            trained = cv.Trained;
            oot = ~cv.ModelParams.Generator.UseObsForIter;
            T = numel(trained);
            cv = []; %#ok<NASGU>
            
            % Check gamma
            if isempty(gamma)
                if isempty(ngamma) || ~isnumeric(ngamma) || ~isscalar(ngamma) || ngamma<=0
                    error(message('stats:ClassificationDiscriminant:cvshrink:BadNGamma'));
                end
                ngamma = ceil(ngamma);
                mingam = zeros(T,1);
                for t=1:T
                    mingam(t) = trained{t}.MinGamma;
                end
                mingam = max(mingam);
                gamma = linspace(mingam,1,ngamma+1);
            elseif ~isnumeric(gamma) || ~isvector(gamma) || any(gamma<0) || any(gamma>1)
                error(message('stats:ClassificationDiscriminant:cvshrink:BadGamma'));
            end
            gamma = sort(gamma);
            Ngamma = numel(gamma);
            gamma = gamma(:);

            % Check delta
            if isempty(delta) || ~isempty(ndelta)
                if isempty(ndelta) || ~isnumeric(ndelta) || ~isscalar(ndelta) || ndelta<0
                    error(message('stats:ClassificationDiscriminant:cvshrink:BadNDelta'));
                end
                ndelta = ceil(ndelta);
                
                % Get delta spacings for every value of gamma
                delta = zeros(Ngamma,ndelta+1);
                for ngam=1:Ngamma
                    maxDelta = max(deltaRange(this.Impl,gamma(ngam)));
                    delta(ngam,:) = linspace(0,maxDelta,ndelta+1);
                end
                
                % Add something small to max delta to make sure that the
                % number of predictors be set to zero at max delta.
                delta(:,end) = delta(:,end) + 10*eps(delta(:,end));
            else
                if ~isnumeric(delta) || ~ismatrix(delta) || any(delta(:)<0)
                    error(message('stats:ClassificationDiscriminant:cvshrink:BadDelta'));
                end
                
                % Interpret row-vector delta as uniform mesh
                if size(delta,1)==1
                    delta = repmat(delta,Ngamma,1);
                else
                    if size(delta,1)~=Ngamma
                        error(message('stats:ClassificationDiscriminant:cvshrink:DeltaSizeMismatch', Ngamma));
                    end
                end
            end

            % Initialize
            Ndelta = size(delta,2);
            K = numel(this.ClassSummary.ClassNames);
            [N,p] = size(this.X);
            mcr = NaN(Ngamma,Ndelta);
            
            % Output number of predictors in the model?
            dopred = nargout>3 && any(delta(:)>0) ...
                && ~isempty(strfind(lower(this.Impl.Type),'linear'));
            npred = repmat(p,Ngamma,Ndelta);
            
            % Loop over gamma and delta
            for ngam=1:Ngamma
                if verbose>0
                    fprintf(1,'%s\n',getString(message('stats:ClassificationDiscriminant:cvshrink_ProcessingGamma',ngam,Ngamma)));
                end
                
                % Set gamma for cross-validated folds. If need to count the
                % number of predictors in the model, set gamma for the
                % object itself.
                if dopred
                    this = setGamma(this,gamma(ngam));
                end
                for t=1:T
                    trained{t} = setGamma(trained{t},gamma(ngam));
                end
                                
                % Loop over CV folds and Delta values. The outer loop is
                % over folds and the inner loop is over Deltas because
                % changing Delta is very quick once the inverse of the
                % covariance matrix has been computed. To avoid
                % re-computing the inverse of the covariance matrix on the
                % fly every time PREDICT is called, turn off memory
                % optimization. Set SaveMemory back to what is used by the
                % CV folds when we no longer need to call PREDICT.
                Sfit = NaN(N,K,Ndelta);
                for t=1:T
                    if verbose>1
                        fprintf(1,'\t\t %s\n',getString(message('stats:ClassificationDiscriminant:cvshrink_ProcessingFold',t)));
                    end
                    if any(delta(:)>0)
                        oldSaveMemory = trained{t}.Impl.SaveMemory;
                        trained{t}.Impl.SaveMemory = false;
                    end
                    for ndel=1:Ndelta
                        trained{t} = setDelta(trained{t},delta(ngam,ndel),false);
                        [~,Sfit(oot(:,t),:,ndel)] = predict(trained{t},this.X(oot(:,t),:));
                    end
                    if any(delta(:)>0)
                        trained{t}.Impl.SaveMemory = oldSaveMemory;
                    end
                end

                % Loop over Delta to compute the number of predictors in
                % the model (if needed) and cross-validated loss.
                % nLinearCoeffs is vectorized wrt delta as long as delta is
                % a vector.
                if dopred
                    npred(ngam,:) = nLinearCoeffs(this,delta(ngam,:));
                end
                for ndel=1:Ndelta
                    mcr(ngam,ndel) = classreg.learning.loss.mincost(C,Sfit(:,:,ndel),this.W,cost);
                end
            end
        end        
    end
            
    methods(Static,Hidden)
        function this = loadobj(obj)
            if isempty(obj.Impl)
                % Load 11b discriminant
                modelParams = fillDefaultParams(obj.ModelParams,...
                    obj.X,obj.PrivY,obj.W,obj.DataSummary,obj.ClassSummary);
                this = ClassificationDiscriminant(obj.X,obj.PrivY,obj.W,...
                    modelParams,obj.DataSummary,obj.ClassSummary,...
                    obj.PrivScoreTransform);
            else
                % Load post-11b discriminant
                this = obj;
            end
        end
    end
    
end
