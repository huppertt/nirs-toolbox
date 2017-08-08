classdef ClassificationModel < classreg.learning.Predictor
%ClassificationModel Compact classification model.
%   ClassificationModel is the super class for compact classification
%   models.

%   Copyright 2010-2014 The MathWorks, Inc.

    properties(GetAccess=public,SetAccess=protected,Hidden=true)
        ClassSummary = struct('ClassNames',{},'NonzeroProbClasses',{},'Cost',[],'Prior',[]);
    end
    
    properties(GetAccess=public,SetAccess=protected,Hidden=true)
        PrivScoreTransform = [];
        PrivScoreType = [];
    end
    
    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        %CLASSNAMES Names of classes in Y.
        %   The ClassNames property is an array containing the class names for the
        %   response variable Y.
        %
        %   See also classreg.learning.classif.ClassificationModel.
        ClassNames;
    end

    properties(GetAccess=public,SetAccess=public,Dependent=true)
        %PRIOR Prior class probabilities.
        %   The Prior property is a vector with prior probabilities for classes.
        %
        %   See also classreg.learning.classif.ClassificationModel.
        Prior;
        
        %COST Misclassification costs.
        %   The Cost property is a square matrix with misclassification costs.
        %   Cost(I,J) is the cost of misclassifying class ClassNames(I) as class
        %   ClassNames(J).
        %
        %   See also classreg.learning.classif.ClassificationModel.
        Cost;

        %SCORETRANSFORM Transformation applied to predicted classification scores.
        %   The ScoreTransform property is a string describing how raw
        %   classification scores predicted by the model are transformed. You can
        %   assign a function handle or one of the following strings to this
        %   property: 'none', 'doublelogit', 'identity', 'invlogit', 'ismax',
        %   'logit', 'sign', 'symmetricismax', 'symmetriclogit', and 'symmetric'.
        %   You can use either 'identity' or 'none' for the identity
        %   transformation.
        %
        %   See also classreg.learning.classif.ClassificationModel.
        ScoreTransform;
    end
       
    properties(GetAccess=public,SetAccess=public,Dependent=true,Hidden=true)
        %ScoreType Type of the score returned by this classifier.
        %   The ScoreType property is a string set to one of:
        %       'probability'     - Classifier score is class posterior
        %                           probability.
        %       '01'              - Classifier score ranges from 0 to 1 and does
        %                           not represent posterior probability.
        %       'inf'             - Classifier score ranges from 0 to +Inf or from
        %                           -Inf to +Inf.
        %       'unknown'         - Classifier score does not have a well-defined
        %                           range, likely because the ScoreTransform
        %                           property does not have one of the expected
        %                           values.
        %   Set ScoreType after setting ScoreTransform. When you assign into the
        %   ScoreTransform property, the ScoreType property is reset.
        ScoreType;
    end
        
    properties(GetAccess=public,SetAccess=protected,Dependent=true,Hidden=true)
        %ContinuousLoss Loss function most appropriate for this classifier.
        %   The ContinuousLoss property is a function handle for one of the
        %   functions under classreg.learning.loss.
        ContinuousLoss;
    end
    
    properties(GetAccess=public,SetAccess=public,Hidden=true)
        DefaultLoss = @classreg.learning.loss.mincost;
        LabelPredictor = @classreg.learning.classif.ClassificationModel.minCost;
        DefaultScoreType = 'probability';
    end
    
    methods
        function cnames = get.ClassNames(this)
            cnames = labels(this.ClassSummary.ClassNames);
        end
        
        function cost = get.Cost(this)
            K = length(this.ClassSummary.ClassNames);
            if isempty(this.ClassSummary.Cost)
                cost = ones(K) - eye(K);
            else
                cost = zeros(K);
                [~,pos] = ismember(this.ClassSummary.NonzeroProbClasses,...
                    this.ClassSummary.ClassNames);
                cost(pos,pos) = this.ClassSummary.Cost;
                unmatched = 1:K;
                unmatched(pos) = [];
                cost(:,unmatched) = NaN;
                cost(1:K+1:end) = 0;
            end
        end

        function prior = get.Prior(this)
            K = length(this.ClassSummary.ClassNames);
            prior = zeros(1,K);
            [~,pos] = ismember(this.ClassSummary.NonzeroProbClasses,...
                this.ClassSummary.ClassNames);
            prior(pos) = this.ClassSummary.Prior;
        end
        
        function this = set.Prior(this,prior)
            this = setPrior(this,prior);
        end
        
        function this = set.Cost(this,cost)
            this = setCost(this,cost);
        end
        
        function st = get.ScoreTransform(this)
            st = classreg.learning.internal.convertScoreTransform(...
                this.PrivScoreTransform,'string',[]);
        end
        
        function this = set.ScoreTransform(this,st)
            this.PrivScoreTransform = ...
                classreg.learning.internal.convertScoreTransform(st,...
                'handle',numel(this.ClassSummary.ClassNames));
            this.PrivScoreType = [];
        end
        
        function scoreType = get.ScoreType(this)
            scoreType = getScoreType(this);
        end

        function this = set.ScoreType(this,st)
            this = setScoreType(this,st);
        end
        
        function cl = get.ContinuousLoss(this)
            cl = getContinuousLoss(this);
        end
    end
    
    methods(Access=protected,Abstract=true)
        s = score(this,X,varargin)
    end
    
    methods(Access=protected)
        function this = setPrior(this,~) %#ok<INUSD>
            error(message('stats:classreg:learning:classif:ClassificationModel:setPrior:Noop'));
        end
        
        function this = setCost(this,~) %#ok<INUSD>
            error(message('stats:classreg:learning:classif:ClassificationModel:setCost:Noop'));
        end
        
        function scoreType = getScoreType(this)            
            if     isequal(this.PrivScoreTransform,@classreg.learning.transform.identity)
                scoreType = this.DefaultScoreType;
            elseif ~isempty(this.PrivScoreType)
                scoreType = this.PrivScoreType;
            else
                scoreType = 'unknown';
            end
        end
        
        function this = setScoreType(this,st)
            this.PrivScoreType = classreg.learning.internal.convertScoreType(st);
        end
        
        function cl = getContinuousLoss(this)
            cl = [];
            if strcmp(this.ScoreType,'probability')
                cl = @classreg.learning.loss.quadratic;
            end
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.Predictor(this,s);
            cnames = this.ClassNames;
            if ischar(cnames)
                s.ClassNames = cnames;
            else
                s.ClassNames = cnames';
            end
            s.ScoreTransform = this.ScoreTransform;
        end
                
        function [labels,posterior,cost] = predictEmptyX(this,X)
            D = numel(this.PredictorNames);
            if size(X,2)~=D
                error(message('stats:classreg:learning:classif:ClassificationModel:predictEmptyX:XSizeMismatch', D));
            end
            labels = repmat(this.ClassNames(1,:),0,1);
            K = numel(this.ClassSummary.ClassNames);
            posterior = NaN(0,K);
            cost = NaN(0,K);
        end        
    end

    methods(Hidden)
        function this = ClassificationModel(dataSummary,classSummary,...
                scoreTransform,scoreType)
            this = this@classreg.learning.Predictor(dataSummary);
            this.ClassSummary = classSummary;
            this.PrivScoreTransform = scoreTransform;
            this.PrivScoreType = scoreType;
        end
        
        function this = setPrivatePrior(this,prior)
            if isempty(prior) || strncmpi(prior,'empirical',length(prior))
                error(message('stats:classreg:learning:classif:ClassificationModel:setPrivatePrior:EmpiricalOrEmptyPrior'));
            end
            this.ClassSummary.Prior = ...
                classreg.learning.classif.FullClassificationModel.processPrior(...
                prior,[],this.ClassSummary.ClassNames,...
                this.ClassSummary.NonzeroProbClasses);
        end
        
        function this = setPrivateCost(this,cost)
            this.ClassSummary.Cost = ...
                classreg.learning.classif.FullClassificationModel.processCost(...
                cost,this.ClassSummary.Prior,this.ClassSummary.ClassNames,...
                this.ClassSummary.NonzeroProbClasses);
        end
        
        function [X,C,W,Y,rowData] = prepareDataForLoss(this,X,Y,W,rowData,cleanRows)
            % Check X
            if ~isnumeric(X) || ~ismatrix(X)
                error(message('stats:classreg:learning:classif:ClassificationModel:prepareDataForLoss:BadXType'));
            end
            
            % Convert to class labels
            Y = classreg.learning.internal.ClassLabel(Y);
            
            % Check size
            N = size(X,1);
            if numel(Y)~=N
                error(message('stats:classreg:learning:classif:ClassificationModel:prepareDataForLoss:SizeXYMismatch'));
            end
            
            % Check weights
            if ~isfloat(W) || ~isvector(W) || length(W)~=N || any(W<0)
                error(message('stats:classreg:learning:classif:ClassificationModel:prepareDataForLoss:BadWeights', N));
            end
            internal.stats.checkSupportedNumeric('Weights',W,true);
            W = W(:);
            
            % Clean rows for missing labels, zero prior and cost?
            if nargin<6 && N>0
                cleanRows = true;
            end
            
            % Any rowData present?
            if nargin>4 && ~isempty(rowData)
                haveRowData = true;
                if size(rowData,1)~=N
                    error(message('stats:classreg:learning:classif:ClassificationModel:prepareDataForLoss:SizeRowDataMismatch',N));
                end
            else
                haveRowData = false;
            end
            
            % Check for missing class labels
            t = ismissing(Y);
            if any(t) && cleanRows
                Y(t) = [];
                X(t,:) = [];
                W(t,:) = [];
                if haveRowData
                    rowData(t,:) = [];
                end
            end
            
            % Get a matrix of class counts.
            % Must use ClassSummary.ClassNames instead of NonzeroProbClasses
            % because the customer can pass a full list of classes.
            C = classreg.learning.internal.classCount(this.ClassSummary.ClassNames,Y);
            
            % Remove observations for classes with zero probability. For
            % one-class learning, Prior is a scalar, and no data reduction
            % should be carried out.
            zeroprior = this.Prior==0;
            if any(zeroprior) && cleanRows && ~isscalar(this.Prior)
                t = any(C(:,zeroprior),2);
                Y(t) = [];
                X(t,:) = [];
                C(t,:) = [];
                W(t,:) = [];
                if haveRowData
                    rowData(t,:) = [];
                end
            end
            
            % Remove observations for classes with zero cost.
            % For one-class learning, Cost is scalar 0, and no data
            % reduction should be carried out.
            zerocost = all(this.Cost==0,2)';
            if any(zerocost) && cleanRows && ~isscalar(this.Cost)
                t = any(C(:,zerocost),2);
                Y(t) = [];
                X(t,:) = [];
                C(t,:) = [];
                W(t,:) = [];
                if haveRowData
                    rowData(t,:) = [];
                end
            end
            
            % Normalize weights to the class prior
            if ~isempty(C)
                WC = bsxfun(@times,C,W);
                Wj = sum(WC,1);
                adjWFactor = zeros(1,numel(Wj));
                zeroprior = Wj==0;
                adjWFactor(~zeroprior) = this.Prior(~zeroprior)./Wj(~zeroprior);
                W = sum(bsxfun(@times,WC,adjWFactor),2);
            end
        end
    end
       
    methods
        function [labels,scores,cost] = predict(this,X,varargin)
        %PREDICT Predict response of the model.
        %   [LABEL,POSTERIOR,COST]=PREDICT(OBJ,X) returns predicted class labels
        %   LABEL, posterior probabilities POSTERIOR and misclassification costs
        %   COST for model OBJ and a matrix of predictors X. X must be a numeric
        %   matrix of size N-by-P, where P is the number of predictors used for
        %   training this model. Classification labels LABEL have the same type as
        %   Y used for training. Posterior probabilities POSTERIOR are an N-by-K
        %   numeric matrix for N observations and K classes. COST is an N-by-K
        %   matrix with predicted misclassification costs per class. The predicted
        %   label is assigned to the class with the minimal misclassification cost.
        %
        %   See also classreg.learning.classif.ClassificationModel.
  
            % Empty data
            if isempty(X)
                [labels,scores,cost] = predictEmptyX(this,X);
                return;
            end
            
            % Get scores from the compact class
            scores = score(this,X,varargin{:});
            
            % Transform scores and find the most probable class
            [labels,scores,cost] = this.LabelPredictor(this.ClassNames,...
                this.Prior,this.Cost,scores,this.PrivScoreTransform);
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
        %   See also classreg.learning.classif.ClassificationModel,
        %   classreg.learning.classif.ClassificationModel/predict.
        
            [X,C] = prepareDataForLoss(this,X,Y,ones(size(X,1),1),[],false);
            [~,Sfit] = predict(this,X,varargin{:});
            m = classreg.learning.loss.classifmargin(C,Sfit);
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
        %       'Weights'   - Observation weights, a numeric vector of length
        %                     size(X,1). By default, all observation weights are
        %                     set to 1. If you supply weights, EDGE computes
        %                     weighted classification edge.
        %
        %   See also classreg.learning.classif.ClassificationModel,
        %   classreg.learning.classif.ClassificationModel/margin.
        
            % Get observation weights
            N = size(X,1);
            args = {'weights'};
            defs = {ones(N,1)};
            [W,~,extraArgs] = internal.stats.parseArgs(args,defs,varargin{:});
            
            % Prepare data
            [X,C,W] = prepareDataForLoss(this,X,Y,W);
            
            % Get observation margins for hypothesis H.
            [~,Sfit] = predict(this,X,extraArgs{:});

            % Check all arguments
            classreg.learning.internal.classifCheck(C,Sfit,W,[]);
            
            % Get edge            
            e = classreg.learning.loss.classifedge(C,Sfit,W,[]);
        end
        
        function l = loss(this,X,Y,varargin)
        %LOSS Classification error.
        %   L=LOSS(MODEL,X,Y) returns classification cost for model MODEL computed
        %   using matrix of predictors X and true class labels Y. X must be a
        %   numeric matrix of size N-by-P, where P is the number of predictors used
        %   for training this model. Y must be of the same type as MODEL.ClassNames
        %   and have N elements.
        %
        %   L=LOSS(MODEL,X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'LossFun'          - Function handle for loss, or string
        %                            representing a built-in loss function.
        %                            Available loss functions for classification:
        %                            'binodeviance', 'classiferror', 'exponential',
        %                            'hinge', and 'mincost'. If you pass a function
        %                            handle FUN, LOSS calls it as shown below:
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
        %                            'mincost'
        %       'Weights'          - Vector of observation weights. By default the
        %                            weight of every observation is set to 1. The
        %                            length of this vector must be equal to the
        %                            number of rows in X.
        %
        %   See also classreg.learning.classif.ClassificationModel,
        %   classreg.learning.classif.ClassificationModel/predict.  
            
            % Get loss function and observation weights
            N = size(X,1);
            args = {       'lossfun' 'weights'};
            defs = {this.DefaultLoss ones(N,1)};
            [funloss,W,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Prepare data
            [X,C,W] = prepareDataForLoss(this,X,Y,W);

            % Check loss function
            funloss = classreg.learning.internal.lossCheck(funloss,'classification');

            % Get observation margins for hypothesis H.
            [~,Sfit] = predict(this,X,extraArgs{:});

            % Check all arguments
            classreg.learning.internal.classifCheck(C,Sfit,W,this.Cost);

            % Get loss
            l = funloss(C,Sfit,W,this.Cost);
        end
        
        function [h,p,err1,err2] = compareHoldout(this,other,X1,X2,Y,varargin)
        %COMPAREHOLDOUT Compare accuracies of two models using new data.
        %   H=COMPAREHOLDOUT(C1,C2,X1,X2,Y) performs a test of the null hypothesis: two
        %   classifiers, C1 applied to X1 and C2 applied to X2, have equal accuracy for
        %   predicting true class labels Y.
        %       H = 0 => Do not reject the null hypothesis at the 5% significance level.
        %       H = 1 => Reject the null hypothesis at the 5% significance level.
        %
        %   Pass C1 and C2 as classifier objects returned by one of the following
        %   functions: fitcdiscr, fitcecoc, fitensemble, fitcknn, fitcnb, fitcsvm, and
        %   fitctree. Pass X1 and X2 as numeric matrices with an equal number of rows.
        %   Pass Y as a categorical array, character array, logical vector, numeric
        %   vector, or cell array of strings. If you pass Y as a character array, it
        %   must have one class label per row. Y must have as many elements as there are
        %   rows in X1 and X2.
        %
        %   [H,P,E1,E2]=COMPAREHOLDOUT(C1,C2,X1,X2,Y) also returns the p-value P of the
        %   test, classification error E1 for C1 and classification error E2 for C2. If
        %   you pass the 'Cost' parameter, E1 and E2 are misclassification costs.
        %
        %   [...]=COMPAREHOLDOUT(C1,C2,X1,X2,Y,'PARAM1',val1,'PARAM2',val2,...)
        %   specifies one or more of the following name/value pairs:
        %       'Alpha'        - Confidence level, a positive scalar. Default: 0.05
        %       'Alternative'  - String indicating the alternative hypothesis, one of:
        %                          * 'unequal' - COMPAREHOLDOUT tests H0: "C1 and C2
        %                                        have equal accuracy" against H1: "C1
        %                                        and C2 have unequal accuracy".
        %                          * 'less'    - COMPAREHOLDOUT tests H0: "C1 is at
        %                                        least as accurate as C2" against H1:
        %                                        "C1 is less accurate than C2".
        %                          * 'greater' - COMPAREHOLDOUT tests H0: "C1 is at
        %                                        most as accurate as C2" against H1: "C1
        %                                        is more accurate than C2".
        %                          Default: 'unequal'
        %       'ClassNames'   - Array of class names. Use the data type that exists in
        %                        Y. You can use this argument to order the classes or
        %                        select a subset of classes. Default: The class names
        %                        that exist in Y.
        %       'Cost'         - Square matrix, where COST(I,J) is the cost of
        %                        classifying a point into class J if its true class is
        %                        I. Alternatively, COST can be a structure S having two
        %                        fields: S.ClassificationCosts containing the cost
        %                        matrix C, and S.ClassNames containing the class names
        %                        and defining the ordering of classes used for the rows
        %                        and columns of the cost matrix. For S.ClassNames use
        %                        the data type that exists in Y. If you pass 'Cost' as a
        %                        numeric matrix, the order of rows and columns matches
        %                        the order defined by 'ClassNames'. If you pass 'Cost',
        %                        COMPAREHOLDOUT can perform an asymptotic two-sided test
        %                        only, that is, the 'Alternative' and 'Test' parameters
        %                        must be set to 'unequal' and 'asymptotic',
        %                        respectively. Default: []
        %       'CostTest'     - String indicating the type of the cost-sensitive test,
        %                        one of: 'likelihood' or 'chisquare'. If 'likelihood',
        %                        COMPAREHOLDOUT uses a likelihood ratio test. If
        %                        'chisquare', COMPAREHOLDOUT uses a chi-square test.
        %                        This parameter is ignored unless you pass a cost matrix
        %                        using the 'Cost' parameter. The chi-square test
        %                        requires an Optimization Toolbox license. Type 'doc
        %                        testcholdout' for info about cost-sensitive tests.
        %       'Test'         - String, one of: 'asymptotic', 'exact', or 'midp'. Let
        %                        N01 be the number of observations misclassified by C1
        %                        and correctly classified by C2. Let N10 be the number
        %                        of observations correctly classified by C1 and
        %                        misclassified by C2. If you set 'Test' to
        %                          * 'asymptotic' - COMPAREHOLDOUT can perform several
        %                                           tests:
        %                               o If you do not pass 'Cost', COMPAREHOLDOUT
        %                                 performs the asymptotic McNemar test assuming
        %                                 (N01-N10)/sqrt(N01+N10) has a normal
        %                                 distribution with zero mean and unit variance.
        %                               o If you pass 'Cost', COMPAREHOLDOUT by default
        %                                 performs a likelihood ratio test assuming
        %                                 twice the log of the likelihood ratio has a
        %                                 chi-square distribution with one degree of
        %                                 freedom. If you set 'CostTest' to 'chisquare',
        %                                 COMPAREHOLDOUT performs a chi-square test
        %                                 assuming the test statistic has a chi-square
        %                                 distribution with one degree of freedom.
        %                          * 'exact'      - COMPAREHOLDOUT performs the exact
        %                                           conditional McNemar test assuming
        %                                           that N01 has a binomial distribution
        %                                           with N01+N10 trials and binomial
        %                                           parameter 0.5.
        %                          * 'midp'       - COMPAREHOLDOUT performs the mid-p
        %                                           value McNemar test. This test uses
        %                                           the same distribution assumptions as
        %                                           the exact conditional test does. To
        %                                           compute the p-value, this test
        %                                           corrects the binomial CDF
        %                                           P(X<=x;N01+N10,0.5) by subtracting
        %                                           half of the binomial probability
        %                                           mass function at x,
        %                                           0.5*P(X=x;N01+N10,0.5).
        %                          Default: 'midp' if 'Cost' is not passed and
        %                                   'asymptotic' otherwise
        %
        %   See also: testcholdout, predict.

            if ~isa(other,'classreg.learning.classif.ClassificationModel')
                error(message('stats:classreg:learning:classif:ClassificationModel:compareHoldout:BadOtherType'));
            end
            
            Y = classreg.learning.internal.ClassLabel(Y);
            N = numel(Y);
            
            if ~ismatrix(X1) || size(X1,1)~=N
                error(message('stats:classreg:learning:classif:ClassificationModel:compareHoldout:BadPredictorMatrix','X1',N));
            end
            if ~ismatrix(X2) || size(X2,1)~=N
                error(message('stats:classreg:learning:classif:ClassificationModel:compareHoldout:BadPredictorMatrix','X2',N));
            end
        
            Yhat1 = predict(this,X1);
            Yhat2 = predict(other,X2);
           
            [h,p,err1,err2] = testcholdout(Yhat1,Yhat2,Y,varargin{:});
        end
    end
    
    methods(Static=true,Hidden=true)
        function [labels,scores,cost,classnum] = ...
                maxScore(classnames,Prior,Cost,scores,scoreTransform)
            scores = scoreTransform(scores);
            N = size(scores,1);
            notNaN = ~all(isnan(scores),2);
            [~,cls] = max(Prior);
            labels = repmat(classnames(cls,:),N,1);
            [~,classnum] = max(scores(notNaN,:),[],2);
            labels(notNaN,:) = classnames(classnum,:);
            cost = Cost(:,classnum)';
        end
        
        function [labels,scores,cost,classnum] = ...
                minCost(classnames,Prior,Cost,posterior,scoreTransform)
            cost = posterior*Cost;
            N = size(posterior,1);
            notNaN = ~all(isnan(cost),2);
            [~,cls] = max(Prior);
            labels = repmat(classnames(cls,:),N,1);
            [~,classnum] = min(cost(notNaN,:),[],2);
            labels(notNaN,:) = classnames(classnum,:);
            scores = scoreTransform(posterior);
        end
    end

end
