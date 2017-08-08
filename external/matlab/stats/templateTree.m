function temp = templateTree(varargin)
%TEMPLATETREE Create a decision tree template.
%   T=TEMPLATETREE() returns a tree template suitable for use in the
%   FITENSEMBLE and FITCECOC functions.
%
%   T=TEMPLATETREE('PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
%       'AlgorithmForCategorical' - Algorithm to find the best split on a
%                                   categorical predictor in data with 3 or
%                                   more classes. Set to one of: 'exact',
%                                   'pullleft', 'pca' or 'ovabyclass'. By
%                                   default, ClassificationTree selects the
%                                   optimal subset of algorithms for each
%                                   split using the known number of classes
%                                   and levels of a categorical predictor.
%                                   For regression and classification with
%                                   two classes, decision tree always
%                                   performs the exact search.
%       'MaxNumCategories'        - ClassificationTree splits a
%                                   categorical predictor using the exact
%                                   search algorithm if the predictor has
%                                   at most MaxNumCategories levels in the
%                                   split node. Otherwise
%                                   ClassificationTree finds the best
%                                   categorical split using one of inexact
%                                   algorithms. Pass this parameter as a
%                                   numeric non-negative scalar. Passing a
%                                   small value can lead to loss of
%                                   accuracy and passing a large value can
%                                   lead to long computation time and
%                                   memory overload. For regression and
%                                   classification with two classes,
%                                   decision tree always performs the exact
%                                   search. Default: 10
%       'MaxNumSplits'            - Maximal number of decision splits (or
%                                   branch nodes) per tree. 
%                                   Default for bagging and ECOC: 
%                                     N-1 for data with N observations
%                                   Default for boosting: 1
%       'MergeLeaves'             - When 'on', decision tree merges
%                                   leaves that originate from the same
%                                   parent node, and that give a sum of
%                                   risk values greater or equal to the
%                                   risk associated with the parent node.
%                                   When 'off', decision tree does not
%                                   merge leaves. Defaults: 'off' for
%                                   ensemble fitting and 'on' for ECOC
%       'MinLeafSize'             - Each leaf has at least 'MinLeafSize'
%                                   observations per tree leaf. If you
%                                   supply both 'MinParentSize' and
%                                   'MinLeafSize', decision tree uses the
%                                   setting that gives larger leaves:
%                                   MinParentSize=max(MinParentSize,2*MinLeafSize).
%                                   Default for bagging and boosting: 
%                                     1 for classification and 5 for regression.
%                                   Default for ECOC: 1
%       'MinParentSize'           - Each splitting node in the tree has at
%                                   least 'MinParentSize' observations. If
%                                   you supply both 'MinParentSize' and
%                                   'MinLeafSize', decision tree uses the
%                                   setting that gives larger leaves:
%                                   MinParentSize=max(MinParentSize,2*MinLeafSize).
%                                   Default for bagging and boosting: 
%                                     2 for classification and 10 for regression. 
%                                   Default for ECOC: 10
%       'NumVariablesToSample'    - Number of predictors to select at random
%                                   for each split. Can be a positive integer
%                                   or 'all'; 'all' means use all available
%                                   predictors. Default for bagging: square
%                                   root of the number of predictors for
%                                   classification and one third of
%                                   predictors for regression. Default for
%                                   boosting and ECOC: 'all'.
%       'Prune'                   - When 'on', unpruned trees are grown and
%                                   the optimal sequence of pruned subtrees
%                                   is computed for each tree.  When 'off',
%                                   unpruned trees are grown and the
%                                   optimal sequences are not computed.
%                                   Default for ensemble fitting: 'off'
%                                   Default for ECOC: 'on'
%       'PruneCriterion'          - String with the pruning criterion,
%                                   'error' or 'impurity' for
%                                   classification and 'mse' (mean squared
%                                   error) for regression. Default: 'error'
%                                   for classification and 'mse' for
%                                   regression
%       'QuadraticErrorTolerance' - Tolerance on quadratic error per
%                                   node for regression trees.
%                                   RegressionTree stops splitting nodes
%                                   when quadratic error per node drops
%                                   below TOLER*QED, where QED is the
%                                   quadratic error for the entire data
%                                   computed before the decision tree is
%                                   grown: QED = NORM(Y-YBAR) with YBAR
%                                   estimated as the average of the input
%                                   array Y. Default: 1e-6.
%       'SplitCriterion'          - Criterion for choosing a split. For
%                                   classification, one of: 'gdi' (Gini's
%                                   diversity index), 'twoing' for the
%                                   twoing rule, or 'deviance' for maximum
%                                   deviance reduction (also known as
%                                   cross-entropy). For regression, 'mse'
%                                   (mean squared error). Default: 'gdi'
%                                   for classification and 'mse' for
%                                   regression
%       'Surrogate'               - 'on', 'off', 'all', or a positive
%                                   integer. When 'on', decision tree finds
%                                   10 surrogate splits at each branch
%                                   node. When set to an integer, decision
%                                   tree finds at most the specified number
%                                   of surrogate splits at each branch
%                                   node. When 'all', decision tree finds
%                                   all surrogate splits at each branch
%                                   node. The 'all' setting can use much
%                                   time and memory. Use surrogate splits
%                                   to improve the tree accuracy for data
%                                   with missing values or to compute
%                                   measures of association between
%                                   predictors. Default: 'off'
%
%   See also ClassificationTree, RegressionTree, fitensemble, fitcecoc,
%   fitctree, fitrtree.

%   Copyright 2013-2014 The MathWorks, Inc.

temp = classreg.learning.FitTemplate.make('Tree',varargin{:});
end
