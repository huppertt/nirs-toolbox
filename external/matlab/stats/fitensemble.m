function obj = fitensemble(X,Y,method,nlearn,learners,varargin)
%FITENSEMBLE Fit an ensemble of learners.
%   ENS=FITENSEMBLE(X,Y,METHOD,NLEARN,LEARNERS) fits an ensemble model ENS,
%   which can be used for making predictions on new data. This ensemble
%   model uses a collection of individual learners such as decision trees.
%   These individual learners are grown from one or more templates
%   specified in LEARNERS.
%
%   X must be an N-by-P matrix of predictors with one row per observation
%   and one column per predictor. Y must be an array of N class labels for
%   classification and a vector of floating-point numbers with N elements
%   for regression. For classification Y can be a categorical array,
%   character array, logical vector, numeric vector, or cell array of
%   strings. If Y is a character array, it must have one class label per
%   row.
%
%   ENS is an object of class ClassificationEnsemble for classification and
%   RegressionEnsemble for regression. If you use one of the following five
%   options, ENS is of class ClassificationPartitionedEnsemble for
%   classification and RegressionPartitionedEnsemble for regression:
%   'CrossVal', 'KFold', 'Holdout', 'Leaveout' or 'CVPartition'.
%
%   METHOD must be a string with one of the following values
%   (case-insensitive):
%      for classification with 2 classes:
%          'AdaBoostM1'
%          'LogitBoost'
%          'GentleBoost'
%          'RobustBoost' (requires an Optimization Toolbox license)
%          'LPBoost' (requires an Optimization Toolbox license)
%          'RUSBoost'
%          'TotalBoost' (requires an Optimization Toolbox license)
%          'Bag'
%          'Subspace'
%      for classification with 3 or more classes:
%          'AdaBoostM2'
%          'LPBoost' (requires an Optimization Toolbox license)
%          'RUSBoost'
%          'TotalBoost' (requires an Optimization Toolbox license)
%          'Bag'
%          'Subspace'
%      for regression:
%          'LSBoost'
%          'Bag'
%
%   NLEARN is the number of ensemble learning cycles to be performed. At
%   every training cycle, FITENSEMBLE loops over all learner templates in
%   LEARNERS and trains one weak learner for every template. The number of
%   trained learners in ENS is equal to NLEARN*numel(LEARNERS). Usually, an
%   ensemble with a good predictive power needs between a few hundred and a
%   few thousand weak learners. You do not have to train an ensemble for
%   that many cycles at once. You can start by growing a few dozen
%   learners, inspect the ensemble performance and, if necessary, train
%   more weak learners using RESUME method of the ensemble.
%
%   LEARNERS is a cell array of weak learner templates or a single weak
%   learner template. You must construct every template by calling TEMPLATE
%   method of the appropriate class. For example, call TEMPLATETREE if you
%   want to grow an ensemble of trees. Usually you need to supply only one
%   weak learner template. If you supply one weak learner with default
%   parameters, you can pass LEARNERS in as a string with the name of the
%   weak learner, for example, 'Tree'. Note that the ensemble performance
%   depends on the parameters of the weak learners and you can get poor
%   performance for weak learners with default parameters.
%
%   Use the following learner names and templates:
%           'Discriminant'          templateDiscriminant
%           'KNN'                   templateKNN
%           'Tree'                  templateTree
%
%   If you bag trees, you must pass argument 'type' which can be either
%   'classification' or 'regression', for example:
%
%   ENS=FITENSEMBLE(X,Y,'Bag',100,'Tree','type','classification')
%
%   If you use method 'Subspace' and set NLEARN to
%   'AllPredictorCombinations', FITENSEMBLE constructs
%   NCHOOSEK(size(X,2),NPredToSample) learners for all possible
%   combinations of NPredToSample predictors. You can use only one learner
%   template in this case. For example:
%
%   ENS=FITENSEMBLE(X,Y,'Subspace','AllPredictorCombinations',...
%           'Discriminant','NPredToSample',NPredToSample)
% 
%   ENS=FITENSEMBLE(X,Y,METHOD,NLEARN,LEARNERS,'PARAM1',val1,'PARAM2',val2,...)
%   specifies optional parameter name/value pairs:
%       'CategoricalPredictors' - List of categorical predictors. Pass
%                                 'CategoricalPredictors' as one of:
%                                   * A numeric vector with indices between
%                                     1 and P, where P is the number of
%                                     columns of X.
%                                   * A logical vector of length P, where a
%                                     true entry means that the
%                                     corresponding column of X is a
%                                     categorical variable.
%                                   * 'all', meaning all predictors are
%                                     categorical.
%                                   * A cell array of strings, where each
%                                     element in the array is the name of a
%                                     predictor variable. The names must
%                                     match entries in 'PredictorNames'
%                                     values.
%                                   * A character matrix, where each row of
%                                     the matrix is a name of a predictor
%                                     variable. The names must match
%                                     entries in 'PredictorNames' values.
%                                     Pad the names with extra blanks so
%                                     each row of the character matrix has
%                                     the same length.
%                                 Default: []
%       'CrossVal'              - If 'on', grows a cross-validated ensemble
%                                 with 10 folds. You can use 'KFold',
%                                 'Holdout', 'Leaveout' and 'CVPartition'
%                                 parameters to override this
%                                 cross-validation setting. You can only
%                                 use one of these four options ('KFold',
%                                 'Holdout', 'Leaveout' and 'CVPartition')
%                                 at a time when creating a cross-validated
%                                 ensemble. As an alternative, you can
%                                 cross-validate later using CROSSVAL
%                                 method for ensemble ENS. Default: 'off'
%       'CVPartition'           - A partition created with CVPARTITION to
%                                 use in cross-validated ensemble. 
%       'Holdout'               - Holdout validation uses the specified
%                                 fraction of the data for test, and uses
%                                 the rest of the data for training.
%                                 Specify a numeric scalar between 0 and 1.
%       'KFold'                 - Number of folds to use in cross-validated
%                                 ensemble, a positive integer. Default: 10
%       'Leaveout'              - Use leave-one-out cross-validation by
%                                 setting to 'on'.
%       'NPredToSample'         - Number of predictors to sample at random
%                                 without replacement for method
%                                 'Subspace', an integer between 1 and
%                                 size(X,2). Default: 1
%       'NPrint'                - Print-out frequency, a positive integer
%                                 scalar. By default, this parameter is set
%                                 to 'off' (no print-outs). You can use
%                                 this parameter to keep track of how many
%                                 weak learners have been trained, so far.
%                                 This is useful when you train ensembles
%                                 with many learners on large datasets.
%       'PredictorNames'        - A cell array of names for the predictor
%                                 variables, in the order in which they
%                                 appear in X. Default: {'x1','x2',...}
%       'ResponseName'          - Name of the response variable Y, a
%                                 string. Default: 'Y'
%       'Resample'              - 'on' or 'off'. If 'on', grows an ensemble
%                                 by resampling. By default this parameter
%                                 is set to 'off' for any type of ensemble
%                                 except 'Bag', and boosting is performed
%                                 by reweighting observations at every
%                                 learning iteration. If you set this
%                                 parameter to 'on' for boosting, the
%                                 ensemble is boosted by sampling training
%                                 observations using updated weights as the
%                                 multinomial sampling probabilities.
%       'FResample'             - Fraction of the training set to be
%                                 selected by resampling for every weak
%                                 learner. A numeric scalar between 0 and
%                                 1; 1 by default. This parameter has no
%                                 effect unless you grow an ensemble by
%                                 bagging or set 'resample' to 'on'.
%       'Replace'               - 'on' or 'off', 'on' by default. If 'on',
%                                 FITENSEMBLE samples with replacement; if
%                                 'off', without. This parameter has no
%                                 effect unless you grow an ensemble by
%                                 bagging or set 'resample' to 'on'. If you
%                                 set 'resample' to 'on' and 'replace' to
%                                 'off', FITENSEMBLE samples training
%                                 observations assuming uniform weights and
%                                 boosts by reweighting observations.
%       'Weights'               - Vector of observation weights, one weight
%                                 per observation. For regression,
%                                 FITENSEMBLE normalizes the weights to add
%                                 up to one. For classification,
%                                 FITENSEMBLE normalizes the weights to add
%                                 up to the value of the prior probability
%                                 in the respective class. Default:
%                                 ones(size(X,1),1)
%
%   For classification ensembles you can specify additional optional
%   name/value pairs:
%       'ClassNames'            - Array of class names. Use the data type
%                                 that exists in Y. You can use this
%                                 argument to order the classes or select a
%                                 subset of classes for training. Default:
%                                 The class names that exist in Y.
%       'Cost'                  - Square matrix, where COST(I,J) is the
%                                 cost of classifying a point into class J
%                                 if its true class is I. Alternatively,
%                                 COST can be a structure S having two
%                                 fields: S.ClassificationCosts containing
%                                 the cost matrix C, and S.ClassNames
%                                 containing the class names and defining
%                                 the ordering of classes used for the rows
%                                 and columns of the cost matrix. For
%                                 S.ClassNames use the data type that
%                                 exists in Y. Default: COST(I,J)=1 if
%                                 I~=J, and COST(I,J)=0 if I=J. FITENSEMBLE
%                                 uses the input cost matrix to adjust the
%                                 prior class probabilities. FITENSEMBLE
%                                 passes the adjusted prior probabilities
%                                 and the default cost matrix to its
%                                 learners.
%                              NOTE: For the 'Bag' method, FITENSEMBLE
%                                    generates bootstrap replicas (in-bag
%                                    samples) by oversampling classes with
%                                    large misclassification costs and
%                                    undersampling classes with small
%                                    misclassification costs. Out-of-bag
%                                    (OOB) samples consequently have fewer
%                                    observations from classes with large
%                                    misclassification costs and more
%                                    observations from classes with small
%                                    misclassification costs. For small
%                                    datasets and highly skewed costs, the
%                                    number of OOB observations per class
%                                    can be very low. The OOB loss then has
%                                    large variance and may be hard to
%                                    interpret.
%       'Prior'                 - Prior probabilities for each class.
%                                 Specify as one of:
%                                   * A string:
%                                     - 'empirical' determines class
%                                       probabilities from class
%                                       frequencies in Y
%                                     - 'uniform' sets all class
%                                       probabilities equal
%                                   * A vector (one scalar value for each
%                                     class)
%                                   * A structure S with two fields:
%                                     S.ClassProbs containing a vector of
%                                     class probabilities, and S.ClassNames
%                                     containing the class names and
%                                     defining the ordering of classes used
%                                     for the elements of this vector.
%                                 If you pass numeric values, FITENSEMBLE
%                                 normalizes them to add up to one.
%                                 Default: 'empirical'
%                              NOTE: For the 'Bag' method, FITENSEMBLE
%                                    generates bootstrap replicas (in-bag
%                                    samples) by oversampling classes with
%                                    large prior probabilities and
%                                    undersampling classes with small prior
%                                    probabilities. Out-of-bag (OOB)
%                                    samples consequently have fewer
%                                    observations from classes with large
%                                    prior probabilities and more
%                                    observations from classes with small
%                                    prior probabilities. For small
%                                    datasets and highly skewed prior
%                                    probabilities, the number of OOB
%                                    observations per class can be very
%                                    low. The OOB loss then has large
%                                    variance and may be hard to interpret.
%
%   For AdaBoostM1, AdaBoostM2, LogitBoost, GentleBoost, RUSBoost and
%   LSBoost you can specify additional parameter name/value pairs:
%       'LearnRate'             - Learning rate for shrinkage, a numeric
%                                 scalar between 0 and 1. By default, the
%                                 learning rate is set to 1, and the
%                                 ensemble learns at the maximal possible
%                                 speed. If you set the learning rate to a
%                                 smaller value, the ensemble requires more
%                                 learning iterations but often achieves a
%                                 better accuracy. A popular choice for
%                                 ensemble grown with shrinkage is 0.1.
%
%   For RUSBoost you can specify additional parameter name/value pairs:
%       'RatioToSmallest'       - Either a numeric scalar or vector with K
%                                 elements for K classes. Every element of
%                                 this vector is the sampling proportion
%                                 for this class with respect to the class
%                                 with fewest observations in Y. If you
%                                 pass a scalar, FITENSEMBLE uses this
%                                 sampling proportion for all classes. For
%                                 example, you have class A with 100
%                                 observations and class B with 10
%                                 observations. If you pass [2 1] for
%                                 'RatioToSmallest', every learner in the
%                                 ensemble is trained on 20 observations of
%                                 class A and 10 observations of class B.
%                                 If you pass 2 or [2 2], every learner is
%                                 trained on 20 observations of class A and
%                                 20 observations of class B. If you pass
%                                 'ClassNames', FITENSEMBLE matches
%                                 elements in the array of class names to
%                                 elements in this vector. Default:
%                                 ones(K,1)
%
%   For LPBoost and TotalBoost you can specify additional parameter
%   name/value pairs:
%       'MarginPrecision'       - Margin precision for corrective boosting
%                                 algorithms (LPBoost and TotalBoost), a
%                                 numeric scalar between 0 and 1. This
%                                 parameter affects the number of boosting
%                                 iterations required for conversion. Use a
%                                 small value to grow an ensemble with many
%                                 learners and use a large value to grow an
%                                 ensemble with few learners. Default: 0.01
%
%   For RobustBoost you can specify additional parameter name/value pairs:
%       'RobustErrorGoal'       - Classification error goal for
%                                 RobustBoost, a numeric scalar between 0
%                                 and 1; 0.1 by default. Usually there is
%                                 an optimal range for this parameter for
%                                 your training data. If you set the error
%                                 goal too low or too high, RobustBoost can
%                                 produce a model with poor classification
%                                 accuracy.
%       'RobustMaxMargin'       - Maximal classification margin for
%                                 RobustBoost in the training set, a
%                                 numeric non-negative scalar. RobustBoost
%                                 minimizes the number of observations in
%                                 the training set with classification
%                                 margins below this threshold. Default: 0
%       'RobustMarginSigma'     - Spread of the distribution of
%                                 classification margins over the training
%                                 set for RobustBoost, a numeric positive
%                                 scalar; 0.1 by default.
%
% Example 1: Train and inspect resubstitution error
%    load ionosphere;
%    ada = fitensemble(X,Y,'AdaBoostM1',100,'Tree');
%    plot(resubLoss(ada,'mode','cumulative'));
%    xlabel('Number of decision trees');
%    ylabel('Resubstitution error');
%
% Example 2: Train and estimate generalization error on a holdout sample
%    load ionosphere;
%    ada = fitensemble(X,Y,'AdaBoostM1',100,'Tree','Holdout',0.5);
%    plot(kfoldLoss(ada,'mode','cumulative'));
%    xlabel('Number of decision trees');
%    ylabel('Holdout error');
% 
% See also templateTree, templateKNN, templateDiscriminant,
% classreg.learning.classif.ClassificationEnsemble,
% classreg.learning.regr.RegressionEnsemble
% classreg.learning.classif.ClassificationEnsemble/resume,
% classreg.learning.regr.RegressionEnsemble/resume,
% classreg.learning.partition.ClassificationPartitionedEnsemble,
% classreg.learning.partition.RegressionPartitionedEnsemble.

%   Copyright 2010-2014 The MathWorks, Inc.


narginchk(5,inf);
if ~ischar(method)
    error(message('stats:fitensemble:MethodNameNotChar'));
end
if ~any(strncmpi(method,classreg.learning.ensembleModels(),length(method)))
    error(message('stats:fitensemble:BadMethod', method));
end
temp = classreg.learning.FitTemplate.make(method,'nlearn',nlearn,'learners',learners,varargin{:});
obj = fit(temp,X,Y);
end
