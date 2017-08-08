function temp = templateKNN(varargin)
%TEMPLATEKNN Create a classification KNN template.
%   T=TEMPLATEKNN() returns a KNN classification template suitable for use
%   in the FITENSEMBLE and FITCECOC functions.
%
%   T=TEMPLATEKNN('PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
%       'NumNeighbors'          - A positive integer specifying the number
%                                 of nearest neighbors in X for classifying
%                                 each point when predicting. Default: 1.
%       'NSMethod'              - Nearest neighbors search method.
%                                 Value is either:
%                                   - 'kdtree' uses a kd-tree to find
%                                     nearest neighbors. 'kdtree' is only
%                                     valid when the distance metric
%                                     is one of the following metrics:
%                                       - 'euclidean'
%                                       - 'cityblock'
%                                       - 'minkowski'
%                                       - 'chebyshev'
%                                   - 'exhaustive' uses the exhaustive
%                                     search algorithm. The distance values
%                                     from all the points in X to each
%                                     point in Y are computed to find
%                                     nearest neighbors.
%                                   Default is 'kdtree' when the number of
%                                   columns of X is not greater  than 10, X
%                                   is not sparse, and the distance metric
%                                   is one of the above 4 metrics;
%                                   otherwise, default is 'exhaustive'.
%       'IncludeTies'           - A logical value. Use true to include all
%                                 neighbors whose distance equal the Kth
%                                 smallest distance. Use false to include
%                                 exactly K nearest neighbors.
%       'DistanceWeight'        - A string or a function handle specifying
%                                 the distance weighting function. The
%                                 choices of the string are:
%                                   - 'equal': Each neighbor gets equal
%                                      weight (default).
%                                   - 'inverse': Each neighbor gets weight
%                                      1/d, where d is the distance between
%                                      this neighbor and the point being
%                                      classified.
%                                   - 'squaredinverse': Each neighbor gets
%                                     weight 1/d^2, where d is the distance
%                                     between this neighbor and the point
%                                     being classified.
%                                   - A distance weighting function
%                                     specified using @. A distance
%                                     weighting function must be of the
%                                     form:
%
%                                     function DW = DISTWGT(D)
%
%                                     taking as argument a matrix D and
%                                     returning a matrix of distance weight
%                                     DW. D and DW can only contain
%                                     non-negative numerical values. DW must
%                                     have the same size as D. DW(I,J) is
%                                     the weight computed based on D(I,J).
%        'BreakTies'            - Method of breaking ties if more than one
%                                 class has the same smallest cost. Choices
%                                 are:
%                                   - 'smallest': Assign the point to the
%                                     class with the smallest index. This
%                                     is default.
%                                   - 'nearest': Assign the point to the
%                                     class of its nearest neighbor.
%                                   -  'random': Randomly pick a class
%                                      out the classes with the smallest
%                                      cost.
%        'Distance'              - A string or a function handle specifying
%                                  the distance metric. The value can be
%                                  one of the following:
%                                  'euclidean'   
%                                     - Euclidean distance (default).
%                                  'seuclidean'
%                                     - Standardized Euclidean distance.
%                                       Each coordinate difference between
%                                       X and a query point is divided by
%                                       an element of  vector S. The
%                                       default value of is the standard
%                                       deviation computed from X,
%                                       S=NANSTD(X). To specify another
%                                       value for S, use the 'Scale'
%                                       argument.
%                                  'cityblock'  
%                                     - City Block distance.
%                                  'chebychev'  
%                                     - Chebychev distance (maximum
%                                       coordinate  difference).
%                                  'minkowski'  
%                                     - Minkowski distance. The default
%                                       exponent is 2. To specify a
%                                       different exponent, use the 'P'
%                                       argument.
%                                  'mahalanobis' 
%                                     - Mahalanobis distance, computed
%                                       using a positive definite
%                                       covariance matrix C. The default
%                                       value of C is the sample covariance
%                                       matrix of X, as computed by
%                                       NANCOV(X). To specify another value
%                                       for C, use the 'Cov' argument.
%                                  'cosine' 
%                                     - One minus the cosine of the
%                                       included angle between observations
%                                       (treated as vectors).
%                                   'correlation' 
%                                     - One minus the sample linear
%                                       correlation between observations
%                                       (treated as sequences of values).
%                                   'spearman'   
%                                     - One minus the sample Spearman's
%                                       rank correlation between
%                                       observations (treated as sequences
%                                       of values).
%                                   'hamming'     
%                                     - Hamming distance, percentage of
%                                       coordinates that differ.
%                                   'jaccard'     
%                                     - One minus the Jaccard coefficient,
%                                       the percentage of nonzero
%                                       coordinates that differ.
%                                   function     
%                                     - A distance function specified using
%                                       @ (for example @DISTFUN). A
%                                       distance function must be of the
%                                       form
%
%                                       function D2 = DISTFUN(ZI, ZJ),
%
%                                       taking as arguments a 1-by-N vector
%                                       ZI containing a single row of X or
%                                       Y, an M2-by-N matrix ZJ containing
%                                       multiple rows of X or Y, and
%                                       returning an M2-by-1 vector of
%                                       distances D2, whose Jth element is
%                                       the distance between the
%                                       observations ZI and ZJ(J,:).
%        'Exponent'             - A positive scalar indicating the
%                                   exponent of Minkowski distance. This
%                                   argument is only valid when
%                                  'Distance' is 'minkowski'. Default: 2.
%        'Cov'                  - A positive definite matrix indicating the
%                                 covariance matrix when computing the
%                                 Mahalanobis distance. This argument is
%                                 only valid when 'Distance' is
%                                 'mahalanobis'. Default is NANCOV(X).
%       'Scale'                 - A vector S containing non-negative
%                                 values, with length equal to the number
%                                 of columns in X. Each  coordinate
%                                 difference between X and a query point is
%                                 divided by the corresponding element of
%                                 S. This argument is only valid when
%                                 'Distance' is 'seuclidean'. Default is
%                                 NANSTD(X).
%       'BucketSize'            - The maximum number of data points in the
%                                 leaf node of the kd-tree (default is 50).
%                                 This argument is only meaningful when
%                                 'NSMethod' is 'kdtree'.
%       'Standardize'           - Logical scalar. If true, standardize X by
%                                 centering and dividing columns by their
%                                 standard deviations. Default: false
%
%   See also ClassificationKNN, fitensemble, fitcecoc, fitcknn.

%   Copyright 2013-2014 The MathWorks, Inc.

classreg.learning.FitTemplate.catchType(varargin{:});
temp = classreg.learning.FitTemplate.make('KNN','type','classification',varargin{:});
end
