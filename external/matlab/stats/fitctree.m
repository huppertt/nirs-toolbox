function this = fitctree(X,Y,varargin)
%FITCTREE Fit a classification decision tree.
%   TREE=FITCTREE(X,Y) returns a classification decision tree for
%   predictors X and class labels Y.
%
%   X must be an N-by-P matrix of predictors with one row per observation
%   and one column per predictor. Y must be an array of N class labels. Y
%   can be a categorical array, character array, logical vector, numeric
%   vector, or cell array of strings. If Y is a character array, it must
%   have one class label per row.
%
%   TREE is a classification tree with binary splits. If you use one of the
%   following five options, TREE is of class
%   ClassificationPartitionedModel: 'CrossVal', 'KFold', 'Holdout',
%   'Leaveout' or 'CVPartition'. Otherwise, TREE is of class
%   ClassificationTree.
%
%   TREE=FITCTREE(X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
%       'AlgorithmForCategorical' - Algorithm to find the best split on a
%                                   categorical predictor in data with 3 or
%                                   more classes. Set to one of: 'exact',
%                                   'pullleft', 'pca' or 'ovabyclass'. By
%                                   default, FITCTREE selects the optimal
%                                   subset of algorithms for each split
%                                   using the known number of classes and
%                                   levels of a categorical predictor. For
%                                   two classes, FITCTREE always performs
%                                   the exact search.
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
%                                 I~=J, and COST(I,J)=0 if I=J.
%       'CrossVal'              - If 'on', grows a cross-validated
%                                 tree with 10 folds. You can use 'KFold',
%                                 'Holdout', 'Leaveout' and 'CVPartition'
%                                 parameters to override this
%                                 cross-validation setting. You can only
%                                 use one of these four options ('KFold',
%                                 'Holdout', 'Leaveout' and 'CVPartition')
%                                 at a time when creating a cross-validated
%                                 tree. As an alternative, you can
%                                 cross-validate later using CROSSVAL
%                                 method for tree TREE. Default: 'off'
%       'CVPartition'           - A partition created with CVPARTITION to
%                                 use in cross-validated tree.
%       'Holdout'               - Holdout validation uses the specified
%                                 fraction of the data for test, and uses
%                                 the rest of the data for training.
%                                 Specify a numeric scalar between 0 and 1.
%       'KFold'                 - Number of folds to use in cross-validated
%                                 tree, a positive integer. Default: 10
%       'Leaveout'              - Use leave-one-out cross-validation by
%                                 setting to 'on'.
%       'MaxNumCategories'      - FITCTREE splits a categorical predictor
%                                 using the exact search algorithm if the
%                                 predictor has at most MaxNumCategories
%                                 levels in the split node. Otherwise
%                                 FITCTREE finds the best categorical split
%                                 using one of inexact algorithms. Pass
%                                 this parameter as a numeric non-negative
%                                 scalar. Passing a small value can lead to
%                                 loss of accuracy and passing a large
%                                 value can lead to long computation time
%                                 and memory overload. Default: 10
%       'MaxNumSplits'          - Maximal number of decision splits (or
%                                 branch nodes) per tree. Default:
%                                 size(X,1)-1
%       'MergeLeaves'           - When 'on', FITCTREE merges leaves that
%                                 originate from the same parent node, and
%                                 that give a sum of risk values greater or
%                                 equal to the risk associated with the
%                                 parent node. When 'off', FITCTREE does
%                                 not merge leaves. Default: 'on'
%       'MinLeafSize'           - Each leaf has at least 'MinLeafSize'
%                                 observations per tree leaf. If you supply
%                                 both 'MinParentSize' and 'MinLeafSize',
%                                 FITCTREE uses the setting that gives
%                                 larger leaves:
%                                 MinParentSize=max(MinParentSize,2*MinLeafSize).
%                                 Default: 1
%       'MinParentSize'         - Each splitting node in the tree has at
%                                 least 'MinParentSize' observations. If
%                                 you supply both 'MinParentSize' and
%                                 'MinLeafSize', FITCTREE uses the setting
%                                 that gives larger leaves:
%                                 MinParentSize=max(MinParentSize,2*MinLeafSize).
%                                 Default: 10
%       'NumVariablesToSample'  - Number of predictors to select at random
%                                 for each split. Can be a positive integer
%                                 or 'all'; 'all' means use all available
%                                 predictors. Default: 'all'
%       'PredictorNames'        - A cell array of names for the predictor
%                                 variables, in the order in which they
%                                 appear in X. Default: {'x1','x2',...}
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
%                                 If you pass numeric values, FITCTREE
%                                 normalizes them to add up to one.
%                                 Default: 'empirical'
%       'Prune'                 - When 'on', FITCTREE grows the unpruned
%                                 tree and computes the optimal sequence of
%                                 pruned subtrees. When 'off', FITCTREE
%                                 grows the tree without pruning. Default:
%                                 'on'
%       'PruneCriterion'        - String with the pruning criterion, either
%                                 'error' or 'impurity'. Default: 'error'
%       'ResponseName'          - Name of the response variable Y, a
%                                 string. Default: 'Y'
%       'ScoreTransform'        - Function handle for transforming scores,
%                                 or string representing a built-in
%                                 transformation function. Available
%                                 functions: 'symmetric', 'invlogit',
%                                 'ismax', 'symmetricismax', 'none',
%                                 'logit', 'doublelogit', 'symmetriclogit',
%                                 and 'sign'. Default: 'none'
%       'SplitCriterion'        - Criterion for choosing a split. One of
%                                 'gdi' (Gini's diversity index), 'twoing'
%                                 for the twoing rule, or 'deviance' for
%                                 maximum deviance reduction (also known as
%                                 cross-entropy). Default: 'gdi'
%       'Surrogate'             - 'on', 'off', 'all', or a positive
%                                 integer. When 'on', FITCTREE finds 10
%                                 surrogate splits at each branch node.
%                                 When set to an integer, FITCTREE finds at
%                                 most the specified number of surrogate
%                                 splits at each branch node. When 'all',
%                                 FITCTREE finds all surrogate splits at
%                                 each branch node. The 'all' setting can
%                                 use much time and memory. Use surrogate
%                                 splits to improve the tree accuracy for
%                                 data with missing values or to compute
%                                 measures of association between
%                                 predictors. Default: 'off'
%       'Weights'               - Vector of observation weights, one weight
%                                 per observation. FITCTREE normalizes the
%                                 weights to add up to the value of the
%                                 prior probability in the respective
%                                 class. Default: ones(size(X,1),1)
%
%   Example: Grow a classification tree for Fisher's iris data.
%       load fisheriris
%       t = fitctree(meas,species,'PredictorNames',{'SL' 'SW' 'PL' 'PW'})
%       view(t)
%
%   See also ClassificationTree,
%   classreg.learning.partition.ClassificationPartitionedModel.

%   Copyright 2013-2014 The MathWorks, Inc.

temp = ClassificationTree.template(varargin{:});
this = fit(temp,X,Y);
end
