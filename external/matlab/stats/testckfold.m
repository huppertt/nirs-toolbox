function [h,p,loss1,loss2] = testckfold(c1,c2,X1,X2,Y,varargin)
%TESTCKFOLD Compare accuracies of two classifiers by repeated cross-validation.
%   H=TESTCKFOLD(C1,C2,X1,X2,Y) performs a cross-validation test of the null
%   hypothesis: two classifiers, C1 and C2, have equal accuracy for predictor
%   matrices X1 and X2 and class labels Y. Classifier C1 is applied to predictor
%   matrix X1, and classifier C2 is applied to predictor matrix X2. H indicates
%   the result of the hypothesis test:
%       H = 0 => Do not reject the null hypothesis at the 5% significance level.
%       H = 1 => Reject the null hypothesis at the 5% significance level.
%
%   Pass C1 and C2 as templates returned by one of the following functions:
%   templateDiscriminant, templateECOC, templateEnsemble, templateKNN,
%   templateNaiveBayes, templateSVM, and templateTree. Alternatively, you can
%   pass C1 or C2 as classification models returned by one of the following
%   functions: fitcdiscr, fitcecoc, fitensemble, fitcknn, fitcnb, fitcsvm, and
%   fitctree. Pass X1 and X2 as an N-by-D matrix of predictors with one row per
%   observation and one column per predictor. Y can be a categorical array,
%   character array, logical vector, numeric vector, or cell array of strings.
%   If Y is a character array, it must have one class label per row.
%
%   TESTCKFOLD performs an R-by-K test by dividing X1, X2 and Y in K partitions,
%   training classifiers on each set of K-1 partitions and evaluating their
%   accuracies on the held-out partition. The 1st classifier is trained on
%   partitions of X1, and the 2nd classifier is trained on partitions of X2.
%   TESTCKFOLD repeats this procedure R times to form a Student t or F statistic
%   with an appropriate number of degrees of freedom. See the 'Test' parameter
%   for more detail.
%
%   [H,P,E1,E2]=TESTCKFOLD(C1,C2,X1,X2,Y) also returns the p-value P of the test,
%   and two R-by-K matrices, E1 and E2, holding classification errors for the
%   first and second classifier.
%
%   [...]=TESTCKFOLD(C1,C2,X1,X2,Y,'PARAM1',val1,'PARAM2',val2,...) specifies
%   one or more of the following name/value pairs:
%       'Alpha'                 - Confidence level, a positive scalar. Default:
%                                 0.05
%       'Alternative'           - String indicating the alternative hypothesis,
%                                 one of:
%                                   * 'unequal' - TESTCKFOLD tests H0: "C1 and
%                                                 C2 have equal accuracy"
%                                                 against H1: "C1 and C2 have
%                                                 unequal accuracy".
%                                   * 'less'    - TESTCKFOLD tests H0: "C1 is at
%                                                 least as accurate as C2"
%                                                 against H1: "C1 is less
%                                                 accurate than C2".
%                                   * 'greater' - TESTCKFOLD tests H0: "C1 is at
%                                                 most as accurate as C2"
%                                                 against H1: "C1 is more
%                                                 accurate than C2".
%                                 Default: 'unequal'
%       'X1CategoricalPredictors' - List of categorical predictors in X1.
%                                   Pass 'X1CategoricalPredictors' as one of:
%                                    * A numeric vector with indices between 1
%                                      and P, where P is the number of columns
%                                      of X1.
%                                    * A logical vector of length P, where a
%                                      true entry means that the corresponding
%                                      column of X1 is a categorical variable.
%                                    * 'all', meaning all predictors in X1 are
%                                      categorical.
%                                   Default: []
%       'X2CategoricalPredictors' - List of categorical predictors in X2.
%                                   Pass 'X2CategoricalPredictors' as one of:
%                                    * A numeric vector with indices between 1
%                                      and P, where P is the number of columns
%                                      of X2.
%                                    * A logical vector of length P, where a
%                                      true entry means that the corresponding
%                                      column of X2 is a categorical variable.
%                                    * 'all', meaning all predictors in X2 are
%                                      categorical.
%                                   Default: []
%       'ClassNames'            - Array of class names. Use the data type
%                                 that exists in Y. You can use this argument to
%                                 order the classes or select a subset of
%                                 classes for training. Default: The class names
%                                 that exist in Y.
%       'Cost'                  - Square matrix, where COST(I,J) is the
%                                 cost of classifying a point into class J if
%                                 its true class is I. Alternatively, COST can
%                                 be a structure S having two fields:
%                                 S.ClassificationCosts containing the cost
%                                 matrix C, and S.ClassNames containing the
%                                 class names and defining the ordering of
%                                 classes used for the rows and columns of the
%                                 cost matrix. For S.ClassNames use the data
%                                 type that exists in Y. If you pass 'Cost' as a
%                                 numeric matrix, the order of rows and columns
%                                 matches the order defined by 'ClassNames'.
%                                 Default: COST(I,J)=1 if I~=J, and COST(I,J)=0
%                                 if I=J.
%       'LossFun'               - Function handle for loss, or string
%                                 representing a built-in loss function.
%                                 Available loss functions are: 'binodeviance',
%                                 'classiferror', 'exponential', and 'hinge'. If
%                                 you pass a function handle FUN, loss calls it
%                                 as shown below:
%                                      FUN(C,S,W,COST)
%                                 where C is an N-by-K logical matrix for N rows
%                                 in X and K classes in the ClassNames property,
%                                 S is an N-by-K numeric matrix, W is a numeric
%                                 vector with N elements, and COST is a K-by-K
%                                 numeric matrix. C has one true per row for the
%                                 true class. S is a matrix of classification
%                                 scores for classes with one row per
%                                 observation. W is a vector of observation
%                                 weights. COST is a matrix of misclassification
%                                 costs. If you pass 'LossFun', TESTCKFOLD
%                                 returns values of the specified loss for E1
%                                 and E2. Default: 'classiferror'
%       'Options'               - A struct that contains options specifying
%                                 whether to use parallel computation when
%                                 training binary learners. This argument can be
%                                 created by a call to STATSET. TESTCKFOLD uses
%                                 the following fields:
%                                      'UseParallel'
%                                      'UseSubstreams'
%                                      'Streams'
%                                 For information on these fields see
%                                 PARALLELSTATS.
%                               NOTE: If 'UseParallel' is TRUE and
%                                     'UseSubstreams' is FALSE, then the length
%                                     of 'Streams' must equal the number of
%                                     workers used by TESTCKFOLD. If a parallel
%                                     pool is already open, this will be the
%                                     size of the parallel pool. If a parallel
%                                     pool is not already open, then MATLAB may
%                                     try to open a pool for you (depending on
%                                     your installation and preferences). To
%                                     ensure more predictable results, it is
%                                     best to use the PARPOOL command and
%                                     explicitly create a parallel pool prior to
%                                     invoking TESTCKFOLD with 'UseParallel' set
%                                     to TRUE.
%       'Prior'                 - Prior probabilities for each class.
%                                 Specify as one of:
%                                   * A string:
%                                     - 'empirical' determines class
%                                       probabilities from class frequencies in
%                                       Y
%                                     - 'uniform' sets all class probabilities
%                                       equal
%                                   * A vector (one scalar value for each class)
%                                   * A structure S with two fields:
%                                     S.ClassProbs containing a vector of class
%                                     probabilities, and S.ClassNames containing
%                                     the class names and defining the ordering
%                                     of classes used for the elements of this
%                                     vector.
%                                 If you pass numeric values, classifiers
%                                 normalize them to add up to one. If you pass
%                                 'Prior' as a numeric vector, the order of
%                                 elements matches the order defined by
%                                 'ClassNames'. Default: 'empirical'
%       'Test'                  - String, one of: '5x2t', '5x2F' or '10x10t'. If
%                                 '5x2t', TESTCKFOLD runs the evaluation
%                                 procedure 5 times with 2 partitions and uses a
%                                 Student t statistic with 5 degrees of freedom.
%                                 If '5x2F', TESTCKFOLD runs the evaluation
%                                 procedure 5 times with 2 partitions and uses
%                                 an F statistic with 10 and 5 degrees of
%                                 freedom. If '10x10t', TESTCKFOLD runs the
%                                 evaluation procedure 10 times with 10
%                                 partitions and uses a Student t statistic with
%                                 10 degrees of freedom. Default: '5x2F'
%       'Verbose'               - Verbosity level, a non-negative integer. Set
%                                 above 0 to see diagnostic messages. Default: 0
%       'Weights'               - Vector of observation weights, one weight
%                                 per observation. Classifiers normalize the
%                                 weights to add up to the value of the prior
%                                 probability in the respective class. Default:
%                                 ones(size(X1,1),1)
%
% Example 1: Compare SVM and bagged trees on ionosphere data. Observe that SVM
%            gives a smaller error on average, but this improvement is not
%            statistically significant.
%   load ionosphere;
%   c1 = templateSVM('Standardize',true,'KernelFunction','RBF','KernelScale','auto');
%   c2 = templateEnsemble('Bag',200,'Tree','Type','classification');
%   rng('default'); % set the RNG seed for reproducibility
%   [h,p,err1,err2] = testckfold(c1,c2,X,X,Y,'Verbose',1)
%   mean(err1(:)-err2(:)) % mean error difference between the two models
%
% Example 2: Test if removing the last two predictors from the Fisher iris data
%            affects the accuracy of classification by linear SVM. Use the
%            one-vs-one approach to train a multiclass SVM model. Execute in
%            parallel.
%   load fisheriris;
%   c = templateECOC;
%   [h,p] = testckfold(c,c,meas,meas(:,1:2),species,'Test','10x10t',...
%       'Options',statset('UseParallel',true))
%
%   See also testcholdout, templateDiscriminant, templateECOC, templateEnsemble,
%   templateKNN, templateNaiveBayes, templateSVM, templateTree, parallelstats,
%   statset.

%   Copyright 2014 The MathWorks, Inc.

% C1 must be either FitTemplate or FullClassificationModel
if     isa(c1,'classreg.learning.FitTemplate')
    c1 = fillIfNeeded(c1,'classification');
elseif isa(c1,'classreg.learning.classif.FullClassificationModel')
    c1 = classreg.learning.FitTemplate.makeFromModelParams(c1.ModelParameters);
else
    error(message('stats:testckfold:BadClassifierObjectType','C1'));
end

% C2 must be either FitTemplate or FullClassificationModel
if     isa(c2,'classreg.learning.FitTemplate')
    c2 = fillIfNeeded(c2,'classification');
elseif isa(c2,'classreg.learning.classif.FullClassificationModel')
    c2 = classreg.learning.FitTemplate.makeFromModelParams(c2.ModelParameters);
else
    error(message('stats:testckfold:BadClassifierObjectType','C2'));
end

% Convert Y to ClassLabel.
Y = classreg.learning.internal.ClassLabel(Y);
nonzeroClassNames = levels(Y);

% For X and Y, check the size only. The rest will be checked by classifiers.
N1 = size(X1,1);
if numel(Y)~=N1
    error(message('stats:testckfold:PredictorMatrixSizeMismatch','X1'));
end

N2 = size(X2,1);
if N1~=N2
    error(message('stats:testckfold:PredictorMatrixSizeMismatch','X2'));
end

% Decode input args
args = {'classnames' 'alpha' 'lossfun' 'alternative' 'test' 'verbose' ...
    'x1categoricalpredictors' 'x2categoricalpredictors' 'cost' 'weights' 'options'};
defs = {          ''    0.05        ''     'unequal' '5x2F'         0 ...
                           []                        []     []        []        []};
[userClassNames,alpha,lossfun,alternative,mode,verbose,cat1,cat2,cost,W,paropts,~,extraArgs] = ...
    internal.stats.parseArgs(args,defs,varargin{:});

% Process weights
if isempty(W)
    W = ones(N1,1);
end
if numel(W)~=N1
    error(message('stats:testckfold:WeightSizeMismatch',N1));
end
W = W(:);


if isempty(userClassNames)
    % If the user has not passed any class names, use those found in the array
    % of true class labels.
    userClassNames = nonzeroClassNames;
else
    userClassNames = classreg.learning.internal.ClassLabel(userClassNames);
    
    % If none of the class names passed by the user is found in the existing
    % class names, error.
    missingC = ~ismember(userClassNames,nonzeroClassNames);
    if all(missingC)
        error(message('stats:classreg:learning:classif:FullClassificationModel:prepareData:ClassNamesNotFound'));
    end
    
    % If the user passed a subset of classes found in the data, remove labels
    % for classes not included in that subset.
    missingC = ~ismember(nonzeroClassNames,userClassNames);
    if any(missingC)
        unmatchedY = ismember(Y,nonzeroClassNames(missingC));
        Y(unmatchedY)    = [];
        X1(unmatchedY,:) = [];
        X2(unmatchedY,:) = [];
        W(unmatchedY)    = [];
        nonzeroClassNames(missingC) = [];
    end
end

% Remove entries for useless classes and check the validity of the cost matrix.
if ~isempty(cost)
    cost = classreg.learning.classif.FullClassificationModel.processCost(...
        cost,ones(1,numel(nonzeroClassNames)),userClassNames,nonzeroClassNames);
end

% If the user has not passed a custom loss function, use classification error.
% Classification error is the only choice for classifiers of different types.
if isempty(lossfun)
    lossfun = 'classiferror';
end
lossfun = classreg.learning.internal.lossCheck(lossfun,'classification');

doclasserr = false;
if isequal(lossfun,@classreg.learning.loss.classiferror)
    doclasserr = true;
end
if ~doclasserr && ~strcmp(c1.Method,c2.Method)
    error(message('stats:testckfold:BadLossFun'));
end

% If a cost matrix has been passed, switch from classification error to minimal
% cost
if doclasserr && ~isempty(cost)
    lossfun = @classreg.learning.loss.mincost;
end

% Check alpha
if ~isscalar(alpha) || ~isfloat(alpha) || ~isreal(alpha) || isnan(alpha) ...
        || alpha<=0 || alpha>=1
    error(message('stats:testckfold:BadAlpha'));
end

% Determine the alternative hypothesis.
alternative = validatestring(alternative,{'unequal' 'less' 'greater'},...
    'testckfold','Alternative');

% Determine the test type. Make sure the test type and alternative are
% compatible.
mode = validatestring(mode,{'5x2t' '5x2F' '10x10t'},'testckfold','Test');

if strcmp(mode,'5x2F') && ~strcmp(alternative,'unequal')
    error(message('stats:testckfold:BadAlternativeTestCombo'));
end

% Check verbosity level.
if ~isscalar(verbose) || ~isnumeric(verbose) || ~isreal(verbose) ...
        || verbose<0 || round(verbose)~=verbose
    error(message('stats:testckfold:BadVerbose'));
end

% Process parallel options
[useParallel,RNGscheme] = ...
    internal.stats.parallel.processParallelAndStreamOptions(paropts,true);

% Set R and K
if     ismember(mode,{'5x2t' '5x2F'})
    R = 5;
    K = 2;
else
    R = 10;
    K = 10;
end


    % Function for computing loss values
    function [l1,l2] = loopBody(r,s)
        if isempty(s)
            s = RandStream.getGlobalStream;
        end
        
        if verbose>0
            fprintf('%s\n',getString(message('stats:testckfold:ReportRunProgress',r,R)));
        end
        
        cvp = cvpartition(Y,'kfold',K,s);
        
        l1 = NaN(1,K);
        l2 = NaN(1,K);

        % Loop over cross-validation folds
        for k=1:K
            if verbose>1
                fprintf('    %s\n',getString(message('stats:testckfold:ReportFoldProgress',k,K)));
            end

            % Indices for training and test
            itrain = training(cvp,k);
            itest = test(cvp,k);
            
            % Train the two models
            m1 = fit(c1,X1(itrain,:),Y(itrain),'categoricalpredictors',cat1,...
                'cost',cost,'weights',W(itrain),extraArgs{:});
            m2 = fit(c2,X2(itrain,:),Y(itrain),'categoricalpredictors',cat2,...
                'cost',cost,'weights',W(itrain),extraArgs{:});
            
            % Get observation weights and true labels for the test data
            w = W(itest);
            y = Y(itest);
            
            if doclasserr
                % Compute classification error or misclassification cost based
                % on predicted labels. Here, we deviate from the classifier
                % objects in the classreg framework which compute classification
                % error and misclassification cost using predicted scores, not
                % labels. Using labels works best for comparing classifiers.
                
                % Get predicted labels
                Yhat1 = classreg.learning.internal.ClassLabel(predict(m1,X1(itest,:)));
                Yhat2 = classreg.learning.internal.ClassLabel(predict(m2,X2(itest,:)));
                
                % Get logical matrix C of size N-by-L for N observations and L
                % classes with class memberships
                C = classreg.learning.internal.classCount(nonzeroClassNames,y);
                
                % Get logical matrices for predicted class labels, similar to C.
                % Yhat1 and Yhat2 cannot have elements not found in
                % nonzeroClassNames because the two classifiers have been
                % trained on Y.
                C1 = classreg.learning.internal.classCount(nonzeroClassNames,Yhat1);
                C2 = classreg.learning.internal.classCount(nonzeroClassNames,Yhat2);

                % Record loss values.
                l1(k) = lossfun(C,C1,w,cost);
                l2(k) = lossfun(C,C2,w,cost);
            else
                % If we are not computing classification error, just use the
                % LOSS method of the classification objects to compute loss
                % values based on scores.
                
                l1(k) = loss(m1,X1(itest,:),y,'lossfun',lossfun,'weights',w);
                l2(k) = loss(m2,X2(itest,:),y,'lossfun',lossfun,'weights',w);
            end
        end
    end

% Compute loss values
[loss1,loss2] = ...
    internal.stats.parallel.smartForSliceout(R,@loopBody,useParallel,RNGscheme);

%
% Analyze computed errors.
%

delta = loss1 - loss2;

% If all loss values are equal, the classifiers are equivalent.
if all( abs(delta(:)) < 100*eps(loss1(:)+loss2(:)) )
    p = 1;
    h = false;
    return;
end

%
% Apply the chosen test.
%

switch mode
    case '5x2t'        
        mdelta_r = mean(delta,2);
        s2_r = sum(bsxfun(@minus,delta,mdelta_r).^2,2);
        s2 = sum(s2_r);
        t = delta(1,1)/sqrt(s2/5);
        
        switch alternative
            case 'unequal'
                p = 2*tcdf(-abs(t),5);
            case 'less'
                % delta has a large positive value under H1
                p = tcdf(t,5,'upper');
            case 'greater'
                % delta has a large negative value under H1
                p = tcdf(t,5);
        end
    
    case '5x2F'        
        mdelta_r = mean(delta,2);
        s2_r = sum(bsxfun(@minus,delta,mdelta_r).^2,2);
        s2 = sum(s2_r);
        F = sum(delta(:).^2)/(2*s2);
        
        p = fcdf(F,10,5,'upper'); % computed only for 'unequal' H1
    
    case '10x10t'        
        m = mean(delta(:));
        s2 = var(delta(:));
        t = m/sqrt(s2/(K+1));
        
        p = tcdf(t,K);

        switch alternative
            case 'unequal'
                p = 2*tcdf(-abs(t),K);
            case 'less'
                % delta has a large positive value under H1
                p = tcdf(t,K,'upper');
            case 'greater'
                % delta has a large negative value under H1
                p = tcdf(t,K);
        end
end

h = p<alpha;

end
