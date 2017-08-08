classdef ClassificationNaiveBayes < ...
        classreg.learning.classif.FullClassificationModel & classreg.learning.classif.CompactClassificationNaiveBayes
%ClassificationNaiveBayes Naive Bayes classification model.
%   ClassificationNaiveBayes is a naive Bayes classification model. This model
%   can predict responses for new data. This model also stores data used for
%   training and can compute resubstitution predictions.
%
%   An object of this class cannot be created by calling the constructor.
%   Use FITCNB to create a ClassificationNaiveBayes object by fitting a
%   naive Bayes model to training data.
%
%   This class is derived from CompactClassificationNaiveBayes.
%
%   ClassificationNaiveBayes properties:
%       Y                       - True class labels used to train this model.
%       X                       - Matrix of predictors used to train this model.
%       W                       - Weights of observations used to train this model.
%       ModelParameters         - Input parameters to the fitting function.
%       NumObservations         - Number of observations.
%       PredictorNames          - Names of predictors used for this model.
%       CategoricalPredictors   - Indices of categorical predictors.
%       ResponseName            - Name of the response variable.
%       ClassNames              - Names of classes in Y.
%       Prior                   - Prior class probabilities.
%       Cost                    - Misclassification costs.
%       ScoreTransform          - Transformation applied to predicted classification scores.
%       DistributionNames       - The distributions used to model each predictor.
%       DistributionParameters  - Estimated parameters for the individual distributions.
%       CategoricalLevels       - The levels of predictors that use the 'mvmn' distribution.
%       Kernel                  - The type of smoothing kernel for each kernel distribution.
%       Support                 - The density support for each kernel distribution.
%       Width                   - The width for each kernel distribution.
%
%   ClassificationNaiveBayes methods:
%       compact                 - Compact this model.
%       compareHoldout          - Compare two models using test data.
%       crossval                - Cross-validate this model.
%       edge                    - Classification edge.
%       loss                    - Classification loss.
%       margin                  - Classification margins.
%       predict                 - Predicted response of this model.
%       logP                    - Log probability of data according to this model.
%       resubEdge               - Resubstitution classification edge.
%       resubLoss               - Resubstitution classification loss.
%       resubMargin             - Resubstitution classification margins.
%       resubPredict            - Resubstitution predicted response.
%
%   Example: Train a naive Bayes model on Fisher iris data
%       load fisheriris
%       nb = fitcnb(meas,species)
%
%   See also fitcnb, classreg.learning.classif.CompactClassificationNaiveBayes.
    
%   Copyright 2014 The MathWorks, Inc.

    properties (GetAccess=protected, SetAccess=protected)
        % PrivW stores the original weights so that we don;t lose them if
        % the user temporarily sets a prior component to 0.
        PrivW = [];
    end
        
    methods(Hidden)
        function this = ClassificationNaiveBayes (X,Y,W,modelParams,...
                dataSummary,classSummary,scoreTransform)
            this = this@classreg.learning.classif.FullClassificationModel(...
                X,Y,W,modelParams,dataSummary,classSummary,scoreTransform);
            this = this@classreg.learning.classif.CompactClassificationNaiveBayes(...
                dataSummary,classSummary,scoreTransform,[],[],[],[]);
            
            % Set and expand DistributionNames. (We can't do it in the
            % compact constructor, even though it's called above.)
            this.DistributionNames = this.ModelParameters.DistributionNames;
            if ischar(this.DistributionNames) && ~strcmp('mn',this.DistributionNames)
                this.DistributionNames = repmat({this.DistributionNames},1,this.NumDims);
            end
            
            % Warn if we find any 'mvmn' predictors that weren't declared
            % categorical, and add them to the 'CategoricalPredictors'
            % property.
            mvmns = find(strcmp('mvmn',this.DistributionNames));
            if any(~ismember(mvmns, this.CategoricalPredictors))
                warning(message('stats:ClassificationNaiveBayes:ClassificationNaiveBayes:SomeMvmnNotCat'));
                this.DataSummary.CategoricalPredictors = union(this.CategoricalPredictors(:)', mvmns);
            end
            
            % Store W into PrivW:
            this.PrivW = W;
            
            % Find and store unique Categorical levels for use in training and prediction
            this.CategoricalLevels = findCategoricalLevels(this.DistributionNames, this.X);
            
            % Train the model on this.X, this.Y and this.W, to get DistributionParameters
            if ~all(strcmp(this.DistributionNames,'mvmn')) && ~isfloat(this.X) % no integers allowed
                internal.stats.checkSupportedNumeric('X',this.X);
            end
            if strcmp(this.DistributionNames, 'mn')
                % 'mn' must be treated as a special case because it uses a
                % single distribution to model all predictor dimensions.
                this.DistributionParameters = fitMNDist(this);
            else
                this.DistributionParameters = fitNonMNDists(this);
            end
        end
    end
    
    methods
        function cmp = compact (this)
            %COMPACT Compact naive Bayes model.
            %   CMP=COMPACT(MODEL) returns an object of class
            %   CompactClassificationNaiveBayes holding the structure of
            %   the trained naive Bayes classifier. The compact object does
            %   not contain X and Y used for training.
            %
            %   See also fitcnb, ClassificationNaiveBayes,
            %   classreg.learning.classif.CompactClassificationNaiveBayes.
            cmp = classreg.learning.classif.CompactClassificationNaiveBayes(...
                this.DataSummary, this.ClassSummary, this.PrivScoreTransform, this.PrivScoreType, ...
                this.DistributionNames, this.DistributionParameters, this.CategoricalLevels);
        end
    end
    
    methods(Static, Hidden)
        %% A static fit method
        function this = fit (X,Y,varargin)
            temp = classreg.learning.FitTemplate.make(...
                'NaiveBayes','type','classification',varargin{:});
            this = fit(temp,X,Y);
        end
        
        %% A template method
        function temp = template (varargin)
            temp = classreg.learning.FitTemplate.make(...,
                'NaiveBayes','type','classification',varargin{:});
        end
    end
    
    %% Model-specific methods:
    methods (Access=protected)
        
        % Set the Prior
        function this = setPrior (this, prior)
            this = setPrior@classreg.learning.classif.CompactClassificationNaiveBayes(this,prior);
            % make the weights sum to the prior within classes
            C = classreg.learning.internal.classCount(this.ClassSummary.ClassNames, this.PrivY);
            WC = bsxfun(@times,C,this.PrivW);
            Wj = sum(WC,1);
            gt0 = Wj>0;
            this.W = sum(bsxfun(@times,WC,this.Prior(gt0)./Wj(gt0)),2);
        end
        
        function distParams = fitMNDist (this)
            % Fit a 'mn' to each class. distParams{c,d} is P(d|c).
            classreg.learning.classif.CompactClassificationNaiveBayes.checkMNData(this.X);  % This is in the compact class because it needs to do the check during prediction, too.
            distParams = cell(this.NumClasses, this.NumDims);
            % remove rows containing any nonfinite elements
            badRows = any(~isfinite(this.X),2);
            X = this.X(~badRows,:);
            W = this.W(~badRows);
            Y = this.PrivY(~badRows);
            ClassCount = classreg.learning.internal.classCount(this.ClassSummary.ClassNames, Y);
            for c = 1:this.NumClasses
                % Fit class c
                if this.NonzeroProbClasses(c);
                    % Do weighted add-1 smoothing
                    posteriorWCounts = ones(1,this.NumDims);                 	% init to the prior for add-1 smoothing.
                    cRows = ClassCount(:, c);
                    if any(cRows)
                        Wc = W(cRows)/sum(W(cRows));                            % renormalize W within class.
                        wCounts = Wc' * X(cRows,:) * sum(cRows);
                        posteriorWCounts = wCounts + posteriorWCounts;
                    end
                    distParams(c,:) = num2cell(posteriorWCounts/sum(posteriorWCounts));	% This row is a vector of Multinomial probabilities.
                else
                    distParams(c,:) = {[]};
                end
            end
        end
        
        function distParams = fitNonMNDists (this)
            % Fit a non-'mn' distribution to each {class, dimension}
            % combination. distParams{c,d} will contain the parameters for
            % the distribution modeling P(X_d|c). An error is thrown if
            % there is not enough data for a [c,d] combination.
            
            % Check all (c,d) combinations for sufficient data
            NoDataCombos = false(this.NumClasses, this.NumDims);
            NonzeroProbClasses = this.NonzeroProbClasses;
            for d=1:this.NumDims
                [C,~,~,labels] = crosstab(isnan(this.X(:,d)), this.Y);
                if size(C,1)>1                                              % There are both NaNs and non-NaNs.
                    NoDataCombos(:,d) = (C(1,:)'==0 & NonzeroProbClasses);
                elseif strcmp(labels{1,1}, '1')                            	% All NaNs.
                    NoDataCombos(:,d) = NonzeroProbClasses;
                end
            end
            % Report (c,d) combinations with no data
            if any(NoDataCombos(:))
                [row, col] = find(NoDataCombos);
                fprintf('  %14s   %14s\n','Class Name','Predictor Name');
                for n=1:numel(col)
                    fprintf('  %14s   %14s\n',...
                        char(this.ClassSummary.ClassNames(row(n))), this.PredictorNames{col(n)} );
                end
                error(message('stats:ClassificationNaiveBayes:ClassificationNaiveBayes:NoDataForUniFit'));
            end
            
            ClassCount = classreg.learning.internal.classCount(this.ClassSummary.ClassNames, this.PrivY);
            distParams = cell(this.NumClasses, this.NumDims);
            for c = 1:this.NumClasses
                if this.NonzeroProbClasses(c);
                    for d = 1:this.NumDims
                        % Get data for this class and dimension
                        cRows = ClassCount(:, c);
                        x = this.X(cRows,d);
                        w = this.W(cRows);
                        % Remove NaNs
                        nanRows = isnan(x);
                        x(nanRows) = [];
                        w(nanRows) = [];
                        % error if zero variance and using normal distribution
                        if strcmp(this.DistributionNames{d},'normal') && var(x)==0
                            error(message('stats:ClassificationNaiveBayes:ClassificationNaiveBayes:ZeroVarianceForUniFit', ...
                                char(this.ClassSummary.ClassNames(c)), this.PredictorNames{d}));
                        end
                        % Get fitting params for this {c,d}, if any
                        [Kernel, Support, Width] = getKernelInputParams(this, c, d);
                        % Do the fit
                        distParams{c,d} = fitUnivariateDist(x, w, this.DistributionNames{d}, ...
                            Kernel, Support, Width, this.CategoricalLevels{d});
                    end
                else
                    distParams(c,:) = {[]};
                end
            end
        end
                
        function L = examplesInClass (this, Y, classNum)
            % Returns a logical array indicating which elements of Y are
            % instances of class number classNum.
            L = ismember(Y, this.ClassSummary.ClassNames(classNum));
        end
        
        function [Kernel, Support, Width] = getKernelInputParams (this, c, d)
            % Return kernel parameters for class c and dimension d
            kernelDims = this.KernelDims;
            if isempty(kernelDims) || ~kernelDims(d)
                Kernel = [];
                Support = [];
                Width = [];
            else
                MP = this.ModelParameters;
                % Kernel may be a string or cell
                if ischar(MP.Kernel)
                    Kernel = MP.Kernel;
                else
                    Kernel = MP.Kernel{d};
                end
                % Support may be a numeric vector, string, or cell.
                if isnumeric(MP.Support) || ischar(MP.Support)
                    Support = MP.Support;
                else
                    Support = MP.Support{d};
                end
                % Width may be a numeric scalar, vector or matrix.
                if isscalar(MP.Width)
                    Width = MP.Width;
                elseif size(MP.Width,1)==1
                    Width = MP.Width(d);
                elseif size(MP.Width,2)==1
                    Width = MP.Width(c);
                else
                    Width = MP.Width(c,d);
                end
            end
        end
        
        %% Properties to display
        function s = propsForDisp (this,s)
            s = propsForDisp@classreg.learning.classif.FullClassificationModel(this,s);
            s = propsForDisp@classreg.learning.classif.CompactClassificationNaiveBayes(this,s);
        end

    end
    
end

%% local functions

function categoricalLevels = findCategoricalLevels (DistributionNames, X)
    % categoricalLevels is always a cell array of length NumDims. For each dimension
    % d that has distribution 'mvmn', categoricalLevels{d} is a numeric vector
    % containing the unique values in column d of X, excluding NaNs. Note
    % that X is numeric.
    D = size(X,2);
    categoricalLevels = cell(1,D);
    if any(strcmp('mvmn', DistributionNames))
        for d = 1:D
            if strcmp('mvmn', DistributionNames{d})
                nanRows = isnan(X(:,d));
                categoricalLevels{d} = unique(X(~nanRows,d));
            end
        end
    end
end

function distParams = fitUnivariateDist (x, w, distName, kernel, support, width, categoricalLevels)
    % Fit a non-'mn' univariate distribution to column vector x, weighted
    % by w.  Preconditions: (1) There is at least 1 data point. (2) x and w
    % contain no NaNs. (3) w need not sum to 1.

    w = w/sum(w);
    % Do the fits
    switch distName     % 'kernel','mvmn','normal'
        case 'kernel'
            if isnan(width) 
                width = [];
            end
            distParams = prob.KernelDistribution.fit(x, ...
                'frequency',w*length(w), ...                            % Pass effective number of counts.
                'kernel',kernel, 'support',support, 'width',width);   	% distParams is a prob.KernelDistribution object.
        case 'mvmn'
            % Do weighted add-1 smoothing
            posteriorCounts = ones(length(categoricalLevels),1);            	% init to the prior for add-1 smoothing.
            [~, levelIndices]= ismember(x, categoricalLevels);                 % Get mvmn Level indices for x. (Indices into categoricalLevels).
            if any(levelIndices)
                wCounts = accumarray(levelIndices, w, size(categoricalLevels))*length(levelIndices);   % Get the effective number of counts for each level.
                posteriorCounts = posteriorCounts + wCounts;
            end
            distParams = posteriorCounts/sum(posteriorCounts);      	% distParams is a vector of Multinomial probabilities.
        case 'normal'
%             distParams = prob.NormalDistribution.fit(x, 'frequency',w); % CAN NOT USE. Rejects real weights.
%             [m,s] = normfit(x,0.05,zeros(size(x)),w*length(x));         % CAN NOT USE. Gives different answer. Probably uses weights differently.
            mu = x'*w;
            sigma = sqrt(classreg.learning.internal.wnanvar(x,w,1)); 	% 1=Unbiased std estimate. (The NaN part is unnecessary.)
            distParams = [mu; sigma];                                	% distParams is [mean, std] of the gaussian fit.
    end
end

