function this = fitcnb(X,Y,varargin)
%FITCNB Fit a naive Bayes classifier to data.
%   MODEL=FITCNB(X,Y) returns a naive Bayes model MODEL for predictors X and
%   class labels Y with two or more classes.
%
%   X must be an N-by-P matrix of predictors with one row per observation
%   and one column per predictor. Y must be an array of N class labels. Y
%   can be a categorical array, character array, logical vector, numeric
%   vector, or cell array of strings. If Y is a character array, it must
%   have one class label per row.
%
%   MODEL is a naive Bayes model. If you use one of the following
%   five options, MODEL is of class ClassificationPartitionedModel:
%   'CrossVal', 'KFold', 'Holdout', 'Leaveout' or 'CVPartition'. Otherwise,
%   MODEL is of class ClassificationNaiveBayes.
%
%   MODEL=FITCNB(X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
%
%       'DistributionNames'  A string or a 1-by-P cell vector of strings,
%                    specifying which distributions FITCNB uses to model
%                    the data. If the value is a string, FITCNB models all
%                    the features using one type of distribution. FITCNB
%                    can also model different features using different
%                    types of distributions. If the value is a cell vector,
%                    its Jth element specifies the distribution FITCNB uses
%                    for the Jth feature. The available types of
%                    distributions are:
%           'normal'  Normal (Gaussian) distribution. 
%           'kernel'  Kernel smoothing density estimate.
%           'mvmn'    Multivariate multinomial distribution for discrete
%                    data. FITCNB assumes each individual feature follows
%                    a multinomial model within a class. The parameters for
%                    a feature include the probabilities of all possible
%                    values that the corresponding feature can take.
%           'mn'      Multinomial distribution for classifying count-
%                    based data such as in the bag-of-tokens model. In the
%                    bag-of-tokens model, the value of the Jth feature is
%                    the number of occurrences of the Jth token in this
%                    observation, so it must be a non-negative integer.
%                    When 'mn' is used, FITCNB considers each
%                    observation as multiple trials of a Multinomial
%                    distribution, and considers each occurrence of a token
%                    as one trial. The number of categories(bins) in this
%                    multinomial model is the number of distinct tokens,
%                    i.e., the number of columns of X.
%           Default: 'mvmn' for variables declared as categorical
%                    via the 'CategoricalPredictors' Name/value pair
%                    below, and 'normal' otherwise.
%       'Kernel'      For use with 'kernel' distributions. The type of
%                    kernel smoother to use. It can be a string or a 1-by-P
%                    cell array of strings.  Each string can be 'normal',
%                    'box', 'triangle', or 'epanechnikov'. If Kernel is a
%                    cell array, then the only elements used are those for
%                    which 'DistributionNames' is 'kernel'.
%                   Default: 'normal'
%       'Support'    For use with 'kernel' distributions. The regions where
%                   the density can be applied.  It can
%                    be a string, a two-element vector as shown below, or a
%                    1-by-P cell array of these values:
%           'unbounded'    The density can extend over the whole
%                         real line.
%           'positive'     The density is restricted to positive values.
%           [L,U]          A two-element vector specifying the finite lower
%                         bound L and upper bound U for the support of the
%                         density.
%                    If Support is a cell array, then the only elements
%                    used are those for which 'DistributionNames' is
%                    'kernel'.
%           Default: 'unbounded'
%       'Width'     For use with 'kernel' distributions. The width of the
%                  kernel smoothing window.  The default is to select a
%                  default width automatically for each combination of
%                  feature and class, using a value that is optimal for a
%                  Gaussian distribution. The value can be specified as one
%                  of the following:
%           scalar          Width for all features in all classes. 
%           row vector      1-by-P vector where the Jth element is the
%                           width for the Jth feature in all classes.
%           column vector	K-by-1 vector where the Ith element specifies the
%                           width for all features in the Ith class. K
%                           represents the number of classes.
%           matrix          K-by-P matrix M where M(I,J) specifies the
%                           width for the Jth feature in the Ith class.
%           Default: Select a width automatically for each combination of
%                    feature and class, using a value that is optimal for a
%                    Gaussian distribution.  If Width is specified and
%                    contains NaNs, then a width is automatically computed
%                    for those NaN entries.
%
%       'CategoricalPredictors' - List of categorical predictors. Pass
%                                 'CategoricalPredictors' as one of:
%                                   * A numeric vector with indices between
%                                     1 and P
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
%                                 the class names that exist in Y
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
%                                 exists in Y. 
%                                 Default: COST(I,J)=1 if
%                                 I~=J, and COST(I,J)=0 if I=J. 
%       'CrossVal'              - If 'on', fits a cross-validated
%                                 model with 10 folds. You can use 'KFold',
%                                 'Holdout', 'Leaveout' and 'CVPartition'
%                                 parameters to override this
%                                 cross-validation setting. You can only
%                                 use one of these four options ('KFold',
%                                 'Holdout', 'Leaveout' and 'CVPartition')
%                                 at a time when creating a cross-validated
%                                 model. As an alternative, you can
%                                 cross-validate later using the CROSSVAL
%                                 method for the naive Bayes model.
%                                 Default: 'off'
%       'CVPartition'           - A partition created with CVPARTITION to
%                                 use in cross-validated naive Bayes. 
%       'Holdout'               - Holdout validation uses the specified
%                                 fraction of the data for test, and uses
%                                 the rest of the data for training.
%                                 Specify a numeric scalar between 0 and 1.
%       'KFold'                 - Number of folds to use in cross-validated
%                                 naive Bayes, a positive integer. Default: 10
%       'Leaveout'              - Use leave-one-out cross-validation by
%                                 setting to 'on'.
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
%                                 If you pass numeric values, FITCNB
%                                 normalizes them to add up to one.
%                                 Default: 'empirical'
%       'ResponseName'          - Name of the response variable Y, a
%                                 string. Default: 'Y'
%       'ScoreTransform'        - Function handle for transforming scores,
%                                 or string representing a built-in
%                                 transformation function. Available
%                                 functions: 'symmetric', 'invlogit',
%                                 'ismax', 'symmetricismax', 'none',
%                                 'logit', 'doublelogit', 'symmetriclogit',
%                                 and 'sign'. Default: 'none'
%       'Weights'               - Vector of observation weights, one weight
%                                 per observation. FITCNB normalizes the
%                                 weights to add up to the value of the
%                                 prior probability in the respective
%                                 class. Default: ones(size(X,1),1)
%
%  FITCNB treats NaNs, empty strings or 'undefined' values as missing
%  values. For missing values in Y, FITCNB removes the corresponding rows
%  of X. For missing values in X, when distribution 'mn' is used, FITCNB
%  removes the corresponding rows of X, otherwise, FITCNB only removes the
%  missing values and uses the values of other features in the
%  corresponding rows of X.
%
%   Example 1: Train a naive Bayes model on data with four classes. Estimate its
%              error by cross-validation.
%       load fisheriris;
%       nb = fitcnb(meas, species);
%       nb = crossval(nb);
%       kfoldLoss(nb)
%
%   Example 2: Train a naive Bayes model on data with four classes, using
%              Normal distributions for predictors 1 and 2 and
%              kernel smoothing distributions for predictors 3 and 4.
%              Estimate its error by cross-validation.
%       load fisheriris;
%       nb = fitcnb(meas, species, 'DistributionNames', {'normal', 'normal', 'kernel', 'kernel'});
%       nb = crossval(nb);
%       kfoldLoss(nb)
%
%   See also ClassificationNaiveBayes,
%   classreg.learning.partition.ClassificationPartitionedModel. 

%   Copyright 2014 The MathWorks, Inc.

this = ClassificationNaiveBayes.fit(X,Y,varargin{:});
end
