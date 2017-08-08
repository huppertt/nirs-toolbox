function this = fitcdiscr(X,Y,varargin)
%FITCDISCR Fit discriminant analysis.
%   DISCR=FITCDISCR(X,Y) returns a discriminant analysis model for
%   predictors X and class labels Y.
%
%   X must be an N-by-P matrix of predictors with one row per observation
%   and one column per predictor. Y must be an array of N class labels. Y
%   can be a categorical array, character array, logical vector, numeric
%   vector, or cell array of strings. If Y is a character array, it must
%   have one class label per row.
%
%   DISCR is a discriminant analysis model. If you use one of the following
%   five options, DISCR is of class ClassificationPartitionedModel:
%   'CrossVal', 'KFold', 'Holdout', 'Leaveout' or 'CVPartition'. Otherwise,
%   DISCR is of class ClassificationDiscriminant.
%
%   DISCR=FITCDISCR(X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
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
%                                 exists in Y. As an alternative, you can
%                                 assign to the Cost property of
%                                 discriminant DISCR. Default: COST(I,J)=1
%                                 if I~=J, and COST(I,J)=0 if I=J.
%       'CrossVal'              - If 'on', creates a cross-validated
%                                 discriminant with 10 folds. You can use
%                                 'KFold', 'Holdout', 'Leaveout' and
%                                 'CVPartition' parameters to override this
%                                 cross-validation setting. You can only
%                                 use one of these four options ('KFold',
%                                 'Holdout', 'Leaveout' and 'CVPartition')
%                                 at a time when creating a cross-validated
%                                 model. As an alternative, you can
%                                 cross-validate later using CROSSVAL
%                                 method for discriminant DISCR. Default:
%                                 'off'
%       'CVPartition'           - A partition created with CVPARTITION to
%                                 use in cross-validated discriminant.
%       'Holdout'               - Holdout validation uses the specified
%                                 fraction of the data for test, and uses
%                                 the rest of the data for training.
%                                 Specify a numeric scalar between 0 and 1.
%       'KFold'                 - Number of folds to use in cross-validated
%                                 discriminant, a positive integer.
%                                 Default: 10
%       'Leaveout'              - Use leave-one-out cross-validation by
%                                 setting to 'on'.
%       'DiscrimType'           - A case-insensitive string with the type
%                                 of the discriminant analysis. Specify as
%                                 one of 'Linear', 'PseudoLinear',
%                                 'DiagLinear', 'Quadratic',
%                                 'PseudoQuadratic' or 'DiagQuadratic'.
%                                 Default: 'Linear'
%       'Gamma'                 - Parameter for regularizing the
%                                 correlation matrix of predictors.
%                                    For linear discriminant, you can pass
%                                    a scalar greater or equal to 0 and
%                                    less or equal to 1. If you pass 0 for
%                                    'gamma' and 'Linear' for
%                                    'DiscrimType', and if the correlation
%                                    matrix is singular, FITCDISCR sets
%                                    'Gamma' to the minimal value required
%                                    for inverting the covariance matrix.
%                                    If you set 'Gamma' to 1, FITCDISCR
%                                    sets the discriminant type to
%                                    'DiagLinear'. If you pass a value
%                                    between 0 and 1, FITCDISCR sets the
%                                    discriminant type to 'Linear'.
%
%                                    For quadratic discriminant, you can
%                                    pass either 0 or 1. If you pass 0 for
%                                    'Gamma' and 'Quadratic' for
%                                    'DiscrimType', and if one of the
%                                    classes has a singular covariance
%                                    matrix, FITCDISCR errors. If you set
%                                    'Gamma' to 1, FITCDISCR sets the
%                                    discriminant type to 'DiagQuadratic'.
%                                 Default: 0
%       'Delta'                 - Threshold on linear coefficients, a
%                                 non-negative scalar. For quadratic
%                                 discriminant, this parameter must be set
%                                 to 0. Default: 0
%       'FillCoeffs'            - If 'on', fills the struct holding
%                                 discriminant coefficients. If 'off', the
%                                 Coeffs property of discriminant DISCR
%                                 will be empty. Default: 'on'
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
%                                 If you pass numeric values, FITCDISCR
%                                 normalizes them to add up to one. As an
%                                 alternative, you can assign to the Prior
%                                 property of discriminant DISCR. Default:
%                                 'empirical'
%       'ResponseName'          - Name of the response variable Y, a
%                                 string. Default: 'Y'
%       'SaveMemory'            - If 'on', defers computing the full
%                                 covariance matrix until it is needed for
%                                 prediction. If 'off', computes and stores
%                                 the full covariance matrix in the
%                                 returned object. Set this parameter to
%                                 'on' if X has thousands of predictors.
%                                 Default: 'off'
%       'ScoreTransform'        - Function handle for transforming scores,
%                                 or string representing a built-in
%                                 transformation function. Available
%                                 functions: 'symmetric', 'invlogit',
%                                 'ismax', 'symmetricismax', 'none',
%                                 'logit', 'doublelogit', 'symmetriclogit',
%                                 and 'sign'. Default: 'none'
%       'Weights'               - Vector of observation weights, one weight
%                                 per observation. Default:
%                                 ones(size(X,1),1)
%
%   See also ClassificationDiscriminant,
%   classreg.learning.partition.ClassificationPartitionedModel.

%   Copyright 2013-2014 The MathWorks, Inc.

this = ClassificationDiscriminant.fit(X,Y,varargin{:});
end
