function temp = templateNaiveBayes(varargin)
%TEMPLATENAIVEBAYES Create a naive Bayes template.
%   T=templateNaiveBayes() returns a naive Bayes template suitable for use
%   in the FITCECOC function.
%
%   T=templateNaiveBayes('PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
%
%      'DistributionNames'  A string or a 1-by-P cell vector of strings,
%                    specifying which distributions are used to model the
%                    data. If the value is a string, all features are
%                    modeled using one type of distribution. Different
%                    features can also be modeled using different types of
%                    distributions. If the value is a cell vector, its Jth
%                    element specifies the distribution used for the Jth
%                    feature. The available types of distributions are:
%          'normal'  Normal (Gaussian) distribution. 
%          'kernel'  Kernel smoothing density estimate.
%          'mvmn'    Multivariate multinomial distribution for discrete
%                    data. Each individual feature is assumed to follow
%                    a multinomial model within a class. The parameters for
%                    a feature include the probabilities of all possible
%                    values that the corresponding feature can take.
%          'mn'      Multinomial distribution for classifying count-
%                    based data such as in the bag-of-tokens model. In the
%                    bag-of-tokens model, the value of the Jth feature is
%                    the number of occurrences of the Jth token in this
%                    observation, so it must be a non-negative integer.
%                    When 'mn' is used, each observation is interpreted as
%                    multiple trials of a Multinomial distribution, and
%                    each occurrence of a token as one trial. The number of
%                    categories(bins) in this multinomial model is the
%                    number of distinct tokens, i.e., the number of columns
%                    of X.
%          Default: 'mvmn' for variables declared as categorical
%                    via the 'CategoricalPredictors' Name/value pair
%                    below, and 'normal' otherwise.
%      'Kernel'      For use with 'kernel' distributions. The type of
%                    kernel smoother to use. It can be a string or a 1-by-P
%                    cell array of strings.  Each string can be 'normal',
%                    'box', 'triangle', or 'epanechnikov'. If Kernel is a
%                    cell array, then the only elements used are those for
%                    which 'DistributionNames' is 'kernel'.
%                   Default: 'normal'
%      'Support'    For use with 'kernel' distributions. The regions where
%                   the density can be applied.  It can
%                    be a string, a two-element vector as shown below, or a
%                    1-by-P cell array of these values:
%          'unbounded'    The density can extend over the whole
%                         real line.
%          'positive'     The density is restricted to positive values.
%          [L,U]          A two-element vector specifying the finite lower
%                         bound L and upper bound U for the support of the
%                         density.
%                    If Support is a cell array, then the only elements
%                    used are those for which 'DistributionNames' is
%                    'kernel'.
%          Default: 'unbounded'
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
%   See also fitcecoc, fitcnb, ClassificationECOC.

%   Copyright 2014 The MathWorks, Inc.

temp = classreg.learning.FitTemplate.make('NaiveBayes',varargin{:});
end
