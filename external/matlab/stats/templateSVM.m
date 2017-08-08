function temp = templateSVM(varargin)
%TEMPLATESVM Create a Support Vector Machine template.
%   T=templateSVM() returns an SVM template suitable for use in the
%   FITCECOC function.
%
%   T=templateSVM('PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
%       'BoxConstraint'         - Positive scalar specifying the box
%                                 constraint. For one-class learning the
%                                 box constraint is always set to 1.
%                                 Default: 1
%       'CacheSize'             - Either positive scalar or string
%                                 'maximal'. If numeric, this parameter
%                                 specifies the cache size in megabytes
%                                 (MB). If set to 'maximal', the cache is
%                                 made large enough to hold the entire Gram
%                                 matrix of size N-by-N for N rows in X.
%                                 Optimizing the cache size can have a
%                                 significant impact on the training speed
%                                 for data with many observations. Default:
%                                 1000
%       'GapTolerance'          - Non-negative scalar specifying tolerance
%                                 for feasibility gap obtained by SMO or
%                                 ISDA. If zero, this parameter is not used
%                                 to check convergence. Default: 0
%       'DeltaGradientTolerance'- Non-negative scalar specifying tolerance
%                                 for gradient difference between upper and
%                                 lower violators obtained by SMO or ISDA.
%                                 If zero, this parameter is not used to
%                                 check convergence. Default: 1e-3 if you
%                                 set 'Solver' to 'SMO' and 0 if you set
%                                 'Solver' to 'ISDA'
%       'KKTTolerance'          - Non-negative scalar specifying tolerance
%                                 for Karush-Kuhn-Tucker (KKT) violation
%                                 obtained by SMO or ISDA. If zero, this
%                                 parameter is not used to check
%                                 convergence. Default: 0 if you set
%                                 'Solver' to 'SMO' and 1e-3 if you set
%                                 'Solver' to 'ISDA'
%       'IterationLimit'        - Positive integer specifying the maximal
%                                 number of iterations for SMO and ISDA.
%                                 SVM training stops when this limit is
%                                 reached, even if optimization did not
%                                 converge. Default: 1e6
%       'KernelFunction'        - String specifying function for computing
%                                 elements of the Gram matrix. Pass as one
%                                 of: 'linear', 'gaussian' (or 'rbf'),
%                                 'polynomial' or name of a function on
%                                 MATLAB path. A kernel function must be of
%                                 the form
%
%                                 function G = KFUN(U, V)
%
%                                 The returned value, G, is a matrix of
%                                 size M-by-N, where M and N are the number
%                                 of rows in U and V, respectively.
%                                 Default: 'linear' for two-class learning
%                                 and 'gaussian' (or 'rbf') for one-class
%                                 learning
%       'KernelScale'           - Either string 'auto' or positive scalar
%                                 specifying the scale factor. If you pass
%                                 'auto', an appropriate scale factor is
%                                 selected using a heuristic procedure. To
%                                 compute the Gram matrix, elements in
%                                 predictor matrix X are divided by this
%                                 factor if the 'KernelFunction' value is
%                                 one of: 'linear', 'gaussian' (or 'rbf'),
%                                 or 'polynomial'. If you pass your own
%                                 kernel function, you must apply scaling
%                                 in that function. Default: 1
%                              NOTE: The heuristic procedure for estimation
%                                    of the scale factor uses subsampling.
%                                    Estimates obtained by this procedure
%                                    can vary from one application to
%                                    another. Set the random number
%                                    generator seed before training for
%                                    reproducibility.
%       'KernelOffset'          - Non-negative scalar. After 
%                                 an element of the Gram matrix is
%                                 computed, this value is added to the
%                                 computed element. Default: 0 if you set
%                                 'Solver' to 'SMO' and 0.1 if you set
%                                 'Solver' to 'ISDA'
%       'PolynomialOrder'       - Positive integer specifying the degree of
%                                 polynomial to be used for polynomial
%                                 kernel. This parameter is used only if
%                                 you set 'KernelFunction' to 'polynomial'.
%                                 Default: 3
%       'NumPrint'              - Non-negative scalar. Diagnostic messages
%                                 are displayed during optimization by SMO
%                                 or ISDA every 'NumPrint' iterations. This
%                                 parameter is used only if you set
%                                 'Verbose' to 1. Default: 1000
%       'OutlierFraction'       - Scalar between 0 (inclusive) and 1
%                                 specifying expected fraction of outlier
%                                 observations in the training set. For
%                                 two-class learning, observations with
%                                 large gradients are removed ensuring that
%                                 the specified fraction of observations be
%                                 removed by the time convergence is
%                                 reached. For one-class learning, the bias
%                                 term is found such that the specified
%                                 fraction of observations in the training
%                                 set has negative scores. Default: 0
%       'Solver'                - String specifying the solver name.
%                                 Specify as one of: 
%                                    * 'SMO'  - Sequential Minimal
%                                               Optimization
%                                    * 'ISDA' - Iterative Single Data
%                                               Algorithm
%                                    * 'L1QP' - L1 soft-margin minimization
%                                               by quadratic programming
%                                               (requires an Optimization
%                                               Toolbox license)
%                                 All solvers implement L1 soft-margin
%                                 minimization. Default: 'ISDA' if you set
%                                 'OutlierFraction' to a positive value for
%                                 two-class learning and 'SMO' otherwise
%       'ShrinkagePeriod'       - Non-negative integer. Observations are
%                                 moved from active set to inactive set
%                                 every 'ShrinkagePeriod' iterations. If
%                                 you pass zero, the active set is not
%                                 shrunk. Shrinking can speed up
%                                 convergence significantly when the
%                                 support vector set is much smaller than
%                                 the training data. If you want to apply
%                                 shrinkage, set 'ShrinkagePeriod' to 1000
%                                 as a rule of thumb. Default: 0
%       'Standardize'           - Logical scalar. If true, standardize X by
%                                 centering and dividing columns by their
%                                 standard deviations. Default: false
%       'SaveSupportVectors'    - Logical scalar. If true, support vectors and
%                                 their Alpha coefficients are saved in a
%                                 compact SVM object. If false, they are not
%                                 saved. You can set this parameter to false for
%                                 linear kernel only. Use this parameter to
%                                 reduce memory consumption by compact linear
%                                 SVM objects. Default: false for ECOC learning
%                                 by linear SVM and true otherwise.
%       'Verbose'               - Verbosity level, one of:
%                                   * 0   - Diagnostic messages are not
%                                           values of convergence criteria
%                                           are not saved.
%                                   * 1   - Diagnostic messages are
%                                           displayed and values of
%                                           convergence criteria are saved
%                                           every 'NumPrint' iterations.
%                                   * 2   - A lot of diagnostic messages is
%                                           displayed and values of
%                                           convergence criteria are saved
%                                           at every iteration.
%
%   See also fitcecoc, fitcsvm, ClassificationECOC, ClassificationSVM.

%   Copyright 2014 The MathWorks, Inc.

temp = classreg.learning.FitTemplate.make('SVM',varargin{:});
end
