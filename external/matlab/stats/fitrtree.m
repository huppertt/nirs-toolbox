function this = fitrtree(X,Y,varargin)
%FITRTREE Fit a regression decision tree.
%   TREE=FITRTREE(X,Y) returns a regression decision tree for predictors X
%   and response Y.
%
%   X must be an N-by-P matrix of predictors with one row per observation
%   and one column per predictor. Y must be a vector of floating-point
%   numbers with N elements.
%
%   FITRTREE grows the tree using MSE (mean squared error) as the splitting
%   criterion.
%
%   TREE is a regression tree with binary splits. If you use one of the
%   following five options, TREE is of class RegressionPartitionedModel:
%   'CrossVal', 'KFold', 'Holdout', 'Leaveout' or 'CVPartition'. Otherwise,
%   TREE is of class RegressionTree.
%
%   TREE=FITRTREE(X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
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
%       'CrossVal'              - If 'on', grows a cross-validated tree
%                                 with 10 folds. You can use 'KFold',
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
%       'MaxNumSplits'          - Maximal number of decision splits (or
%                                 branch nodes) per tree. Default:
%                                 size(X,1)-1
%       'MergeLeaves'           - When 'on', FITRTREE merges
%                                 leaves that originate from the same
%                                 parent node, and that give a sum of risk
%                                 values greater or equal to the risk
%                                 associated with the parent node. When
%                                 'off', FITRTREE does not merge leaves.
%                                 Default: 'on'
%       'MinLeafSize'           - Each leaf has at least 'MinLeafSize'
%                                 observations per tree leaf. If you supply
%                                 both 'MinParentSize' and 'MinLeafSize',
%                                 FITRTREE uses the setting that gives
%                                 larger leaves:
%                                 MinParentSize=max(MinParentSize,2*MinLeafSize).
%                                 Default: 1
%       'MinParentSize'         - Each splitting node in the tree has at
%                                 least 'MinParentSize' observations. If
%                                 you supply both 'MinParentSize' and
%                                 'MinLeafSize', FITRTREE uses the setting
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
%       'Prune'                 - When 'on', FITRTREE grows
%                                 the unpruned tree and computes the
%                                 optimal sequence of pruned subtrees. When
%                                 'off' FITRTREE grows the tree without
%                                 pruning. Default: 'on'
%       'QuadraticErrorTolerance' - Defines tolerance on quadratic error per
%                                   node for regression trees. Splitting
%                                   nodes stops when quadratic error per
%                                   node drops below TOLER*QED, where QED
%                                   is the quadratic error for the entire
%                                   data computed before the decision tree
%                                   is grown: QED = NORM(Y-YBAR) with YBAR
%                                   estimated as the average of the input
%                                   array Y. Default = 1e-6.
%       'ResponseName'          - Name of the response variable Y, a
%                                 string. Default: 'Y'
%       'Surrogate'             - 'on', 'off', 'all', or a positive
%                                 integer. When 'on', FITRTREE finds 10
%                                 surrogate splits at each branch node.
%                                 When set to an integer, FITRTREE finds at
%                                 most the specified number of surrogate
%                                 splits at each branch node. When 'all',
%                                 FITRTREE finds all surrogate splits at
%                                 each branch node. The 'all' setting can
%                                 use much time and memory. Use surrogate
%                                 splits to improve the tree accuracy for
%                                 data with missing values or to compute
%                                 measures of association between
%                                 predictors. Default: 'off'
%       'Weights'               - Vector of observation weights, one weight
%                                 per observation. FITRTREE normalizes the
%                                 weights to add up to one. Default:
%                                 ones(size(X,1),1)
%
%   Example: Grow a regression tree for car data.
%       load carsmall
%       t = fitrtree([Weight Horsepower],MPG,'PredictorNames',{'Weight' 'Horsepower'})
%       view(t)
%
%   See also RegressionTree,
%   classreg.learning.partition.RegressionPartitionedModel.

%   Copyright 2013-2014 The MathWorks, Inc.

temp = RegressionTree.template(varargin{:});
this = fit(temp,X,Y);
end
