classdef CompactClassificationNaiveBayes < classreg.learning.classif.ClassificationModel
%CompactClassificationNaiveBayes Naive Bayes classification model.
%   CompactClassificationNaiveBayes is a naive Bayes classification model. This model
%   can predict responses for new data. 
%
%   CompactClassificationNaiveBayes properties:
%       DistributionNames       - The distributions used to model each predictor.
%       DistributionParameters  - Estimated parameters for the individual distributions.
%       CategoricalLevels       - The levels of predictors that use the 'mvmn' distribution.
%       Kernel                  - The type of smoothing kernel for each kernel distribution.
%       Support                 - The density support for each kernel distribution.
%       Width                   - The width for each kernel distribution.
%       ClassNames              - Names of classes.
%       Prior                   - Prior class probabilities.
%       Cost                    - Misclassification costs.
%       ScoreTransform          - Transformation applied to predicted classification scores.
%       PredictorNames          - Names of predictors used for this model.
%       CategoricalPredictors   - Indices of categorical predictors.
%       ResponseName            - Name of the response variable.
%
%   CompactClassificationNaiveBayes methods:
%       compareHoldout          - Compare two models using test data.
%       edge                    - Classification edge.
%       logP                    - Log probability of data according to this model.
%       loss                    - Classification loss.
%       margin                  - Classification margins.
%       predict                 - Predicted response of this model.
%
%   See also ClassificationNaiveBayes.
    
%   Copyright 2014 The MathWorks, Inc.

    %% Public configuration of the trained model
    properties (GetAccess=public, SetAccess=protected)
        %DISTRIBUTIONNAMESS Distribution names.
        %    The DistributionNames property is a string or a 1-by-P cell array
        %    of strings indicating the types of distributions for all the
        %    predictors. If the model uses the 'mn' distribution, then
        %    DistributionNames is the single string 'mn'. Otherwise
        %    DistributionNames is a cell array and DistributionNames{J} indicates
        %    the distribution type used for the Jth predictor.
        %
        %    The valid single string for this property is:
        %
        %       'mn'     - Multinomial bag-of-tokens model.
        %
        %    The valid strings for in the cell array are:
        %
        %       'normal' - Normal distribution.
        %       'kernel' - Kernel smoothing density estimate.
        %       'mvmn'   - Multivariate multinomial distribution.
        %
        %    See also CLASSIFICATIONNAIVEBAYES.
        DistributionNames = [];
        
        %DISTRIBUTIONPARAMETERS Distribution parameter estimates
        %   The DistributionParameters property is an K-by-P cell array
        %   containing the parameter estimates for the individual distributions.
        %   DistributionParameters{I,J} contains the parameter estimates for the
        %   Jth predictor in the Ith class. DistributionParameters{I,J} is empty if
        %   the Ith class is empty.
        %
        %   The entry in DistributionParameters{I,J} depends on the distribution
        %   type used for the Jth predictor, as follows:
        %     'normal' - A vector of length two. The first element is the
        %                mean, and the second element is standard deviation.
        %     'kernel' - A prob.KernelDistribution object
        %     'mvmn'   - A vector containing the probability for each possible
        %                level of the Jth predictor in the Ith class. The order of
        %                the probabilities is decided by the sorted order of all
        %                the unique levels of the Jth predictor, and is recorded in
        %                the CategoricalLevels property. The probability of the Jth
        %                predictor having level L in class I,
        %
        %                       Prob(predictor J = L | class I)
        %
        %                is estimated as:
        %
        %                   (1 + the weighted number of observations for which
        %                    predictor J = L in class I) /
        %                  (the number of distinct levels in predictor J +
        %                   the weighted number of observations in class I)
        %
        %                The "weighted number of observations" meeting some
        %                condition is defined as the sum of the weights of the
        %                observations meeting the condition, times the size of the
        %                dataset. The sum of the weights over the entire dataset is
        %                automatically normalized to sum to 1.
        %
        %     'mn'     - A scalar representing the probability of the Jth token
        %                appearing in the Ith class, Prob(token J | class I). It is
        %                estimated as:
        %                  (1 + the weighted number of occurrences of token J in
        %                  class I) /  (P + the total weighted number of
        %                  occurrences of all tokens in class I)
        %
        %    See also CLASSIFICATIONNAIVEBAYES.
        DistributionParameters = [];
        
        %CategoricalLevels     Multivariate multinomial levels.
        %   The CategoricalLevels property is a cell array of length P. If
        %   DistributionNames{d} is 'mvmn', then CategoricalLevels{d} is a
        %   numeric vector containing the unique levels of predictor d
        %   found in the training data; otherwise CategoricalLevels{d} is
        %   empty.
        CategoricalLevels = [];
    end
    
    properties (GetAccess=public, SetAccess=protected, Dependent=true)
        %KERNEL     The type of kernel smoother used for each predictor using the
        %   'kernel' distribution. A 1-by-P cell array of strings.  Each string can
        %   be 'normal', 'box', 'triangle', 'epanechnikov', or [] for non-kernel
        %   predictors.
        Kernel
        
        %SUPPORT   The kernel density support regions. (The regions where the
        %   density can be applied).  A 1-by-P cell array of these values:
        %       'unbounded'    (default) The density can extend over the whole real
        %                       line.
        %       'positive'     The density is restricted to positive values.
        %       [L,U]          A two-element vector specifying the lower bound L
        %                      and upper bound U for the support of the density.
        %       []             For non-kernel predictors.
        Support
        
        %WIDTH      The width of the kernel smoothing window.  A K-by-P matrix M
        %   where M(I,J) specifies the width for the Jth predictor in the Ith
        %   class. NaN if predictor J is a non-kernel predictor.
        Width
    end
    
    properties (GetAccess=protected, SetAccess=protected, Dependent=true)
        % A logical vector indicating kernel dimensions.
        KernelDims
        
        % A logical vector indicating which classes have nonzero prob.
        NonzeroProbClasses
    end
    
    properties (Hidden, Dependent=true)
        % These are for convenience
        NumClasses        
        NumDims
        % The following two are aliases to the above used by the Legacy
        % tests, and may be removed when those tests are updated.
        NClasses        
        NDims
    end
    
    methods
        function [labels, posterior, cost] = predict (this, X, varargin)
            %PREDICT Predict response of the model.
            %   [LABEL,POSTERIOR,COST]=PREDICT(NB,X) returns predicted class
            %   labels LABEL, posterior probabilities POSTERIOR and
            %   misclassification costs COST for model NB and a matrix of
            %   predictors X. X must be a numeric matrix of size N-by-P, where
            %   P is the number of predictors used for training this model.
            %   Classification labels LABEL have the same type as Y used for
            %   training. Posterior probabilities POSTERIOR are an N-by-K
            %   numeric matrix for N observations and K classes. COST is an
            %   N-by-K matrix with predicted misclassification costs per class.
            %   The predicted label is assigned to the class with the minimal
            %   misclassification cost.
            
            % Empty data
            if isempty(X)
                [labels,posterior,cost] = predictEmptyX(this,X);
                return;
            end
            
            % Get scores (log P(x,c)), and turn into posteriors P(c|x)
            posterior = softmaxRows(score(this,X,varargin{:}));
            
            % Transform scores and find the most probable class
            [labels,posterior,cost] = this.LabelPredictor(this.ClassNames,...
                this.Prior,this.Cost,posterior,this.PrivScoreTransform);
        end
        
        function L = logP (this, X)
            %LOGP Log of the unconditional probability density.
            %   LOGP = LOGP(NB,X) returns a column vector with N elements
            %   for N rows in matrix of predictors X. Element I in this
            %   vector is the log of the probability of observing row I of
            %   X computed using the naive Bayes model NB. LOGP is NaN for
            %   rows of X containing any NaNs.
            
            % Sum the joint prob over classes, then take the log. Do it
            % stably by subtracting the max score from each row first.
            S = score(this,X);
            M = max(S,[],2);
            L = M + log(sum(exp(bsxfun(@minus, S, M)),2));
            
            % Set components to NaN if X contained a NaN
            L(any(isnan(X),2)) = NaN;
        end
        
        % Dependent property getters
        function N = get.NumClasses (this)
            N = numel(this.ClassSummary.ClassNames);
        end
        
        function N = get.NumDims (this)
            N = numel(this.PredictorNames);
        end
        
        function N = get.NClasses (this)
            % This may be removed when legacy tests are updated.
            N = this.NumClasses;
        end
        
        function N = get.NDims (this)
            % This may be removed when legacy tests are updated.
            N = this.NumDims;
        end
        
        function K = get.Kernel (this)
            % Returns [] if distribution is 'mn'. Otherwise,
            % Returns a 1xNDims cell array: Kernel types obtained from the
            % first nonempty row of the DistributionParameters property.
            if strcmp(this.DistributionNames, 'mn')
                K = [];
            else
                K = cell(1,this.NumDims);
                kernelDims = this.KernelDims;
                row = find(cellfun(@(p) ~isempty(p),this.DistributionParameters(:,1)),1,'first');
                K(kernelDims) = cellfun(@(p) p.Kernel, this.DistributionParameters(row,kernelDims), 'UniformOutput',false);
            end
        end
        
        function S = get.Support (this)
            % Returns [] if distribution is 'mn'. Otherwise,
            % Returns a 1xNDims cell array: Kernel support obtained from
            % the first nonempty row of the DistributionParameters property.
            if strcmp(this.DistributionNames, 'mn')
                S = [];
            else
                S = cell(1,this.NumDims);
                kernelDims = this.KernelDims;
                row = find(cellfun(@(p) ~isempty(p),this.DistributionParameters(:,1)),1,'first');
                S(kernelDims) = cellfun(@(p) p.Support.range, this.DistributionParameters(row,kernelDims), 'UniformOutput',false);
            end
        end
            
        function W = get.Width (this)
            % Returns [] if distribution is 'mn'. Otherwise,
            % Returns a K x P numeric array: Kernel width obtained
            % from the DistributionParameters property.
            if strcmp(this.DistributionNames, 'mn')
                W = [];
            else
                W = NaN(this.NumClasses, this.NumDims);
                nonemptyKernelCells = repmat(this.KernelDims, this.NumClasses, 1) & cellfun(@(p) ~isempty(p),this.DistributionParameters);
                W(nonemptyKernelCells) = cellfun(@(p) double(p.BandWidth), this.DistributionParameters(nonemptyKernelCells));   % double() is required here because the BandWidth is sometimes of type 'single'.
            end
        end
        
        function L = get.KernelDims (this)
            % Return logical vector indicating kernel dimensions.
            L = strcmp('kernel',this.DistributionNames);
        end
        
        function L = get.NonzeroProbClasses (this)
            % Returns a logical vector indicating which classes have nonzero prob.
            L = ismember(this.ClassSummary.ClassNames, this.ClassSummary.NonzeroProbClasses);
        end
    end
    
    methods (Static, Access=protected)
        function checkMNData (X)
            % Error if any finite element of X is negative or non-integral.
            % Rows containing NaNs and Infs will be removed later.
            X = X(:);
            if any(isfinite(X) & (X<0 | floor(X)~=X))
                error(message('stats:ClassificationNaiveBayes:ClassificationNaiveBayes:BadDataForMN'));
            end
        end
    end
    
    methods (Access=protected)
        %% The constructor
        function this = CompactClassificationNaiveBayes (dataSummary,classSummary,scoreTransform,scoreType, ...
                DistributionNames, DistributionParameters, CategoricalLevels)
            this = this@classreg.learning.classif.ClassificationModel(dataSummary,classSummary,scoreTransform,scoreType);
            % Set params
            this.DistributionNames = DistributionNames;
            this.DistributionParameters = DistributionParameters;
            this.CategoricalLevels = CategoricalLevels;
        end
        
     	%% The score method
        function logPxc = score (this, X)
            % Returns an NxC matrix of the log joint probability
            % log(P(x,c)). This can be used later to compute either the
            % posterior P(c|x) or the data likelihood P(x). (But note that
            % if there are missing values, then P(x) is really the marginal
            % probability of the non-missing components of x. This should
            % not be compared to the probability of fully-visible x.)
            
            if ~isfloat(X) || ~ismatrix(X)
                error(message('stats:ClassificationNaiveBayes:ClassificationNaiveBayes:BadX'));
            end
            
            if size(X,2)~=this.NumDims
                error(message('stats:ClassificationNaiveBayes:ClassificationNaiveBayes:BadXSize',this.NumDims));
            end
            
            logPxc = zeros(size(X,1), this.NumClasses);   % preallocate.
            
            % Handle 'mn' special case
            if ischar(this.DistributionNames) && strcmp(this.DistributionNames, 'mn')
                classreg.learning.classif.CompactClassificationNaiveBayes.checkMNData(X);
                for c = 1:this.NumClasses
                    if this.NonzeroProbClasses(c);
                        logPxc(:,c) = classreg.learning.internal.mnlogpdf(X, cell2mat(this.DistributionParameters(c,:))) ...	% log P(X|c)
                                      + log(this.Prior(c));                                         % + log P(c)
                    else
                        % Class has zero prob, so set joint P(x,c)=0.
                        logPxc(:,c) = -Inf;
                    end
                end
            else
                % Handle the non-'mn' cases
                for c = 1:this.NumClasses
                    if this.NonzeroProbClasses(c);
                        logPxdgc = zeros(size(X,1), this.NumDims);    % preallocate.
                        for d = 1:this.NumDims
                            logPxdgc(:,d) = univariateLogP(X(:,d), this.DistributionNames{d}, ...
                                this.DistributionParameters{c,d}, this.CategoricalLevels{d});    % Fill in column d.
                        end
                        % To get log joint, sum non-NaN log probs over d, then add log prior
                        logPxc(:,c) = nansum(logPxdgc, 2) + log(this.Prior(c));   
                    else
                        % Class has zero prob, so set joint P(x,c)=0.
                        logPxc(:,c) = -Inf;
                    end
                end
            end
        end
        
        %% Making the Prior and Cost properties settable
        function this = setPrior (this, prior)
            this = setPrivatePrior(this,prior);
        end
            
        function this = setCost (this, cost)
            this = setPrivateCost(this,cost);
        end
            
        %% Properties to display
        function s = propsForDisp (this,s)
            s = propsForDisp@classreg.learning.classif.ClassificationModel(this,s);
            s.DistributionNames = this.DistributionNames;
            s.DistributionParameters = this.DistributionParameters;
            if any(strcmp('mvmn',this.DistributionNames))
                s.CategoricalLevels = this.CategoricalLevels;
            end
            if any(this.KernelDims)
                s.Kernel = this.Kernel;
                s.Support = this.Support;
                s.Width = this.Width;
            end
        end
    end
end

%% local functions

function Y = softmaxRows (X)
    % Exponentiate then normalize each row of X to sum to 1
    numerators = exp(bsxfun(@minus, X, max(X,[],2)));       % Subtract the max of each row before exp for numerical stability.
    Y = bsxfun(@rdivide, numerators, sum(numerators,2));
end

function logPx = univariateLogP (x, DistName, DistParams, categoricalLevels)
    % Return vector of the log prob of the elements of vector x according
    % to distribution 'Distname', with parameters 'DistParams'.
    switch DistName
        case 'kernel'
            logPx = log(DistParams.pdf(x));                             % DistParams is a prob.KernelDistribution object.
        case 'mvmn'
            logPx = -Inf(size(x));                                      % logP will be -Inf in a row if the value is not found in categoricalLevels...
            logPx(isnan(x)) = NaN;                                      % ... unless the value is NaN, in which case logPx is NaN.
            [found, levelIDs] = ismember(x, categoricalLevels);         % Get mvmn Level IDs from x. (Indices into categoricalLevels).
            logPx(found) = log(DistParams(levelIDs(found)));            % For every levelID, get the log of its DistParam.
        case 'normal'
            logPx = classreg.learning.internal.normlogpdf(x, DistParams(1), DistParams(2));
    end
end
