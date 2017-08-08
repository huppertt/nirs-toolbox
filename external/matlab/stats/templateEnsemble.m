function temp = templateEnsemble(method,nlearn,learners,varargin)
%TEMPLATEENSEMBLE Create a template for ensemble learning.
%   T=templateEnsemble(METHOD,NLEARN,LEARNERS) returns an ensemble template
%   suitable for use in the FITCECOC function.
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
%   every training cycle, a loop is carried over all learner templates in
%   LEARNERS and one weak learner is trained for every template. The number
%   of trained learners in ENS is equal to NLEARN*numel(LEARNERS). Usually,
%   an ensemble with a good predictive power needs between a few hundred
%   and a few thousand weak learners. 
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
%   If you use method 'Subspace' and set NLEARN to
%   'AllPredictorCombinations', NCHOOSEK(size(X,2),NPredToSample) learners
%   are constructed for all possible combinations of NPredToSample
%   predictors. You can use only one learner template in this case. For
%   example:
%
%   T=templateEnsemble('Subspace','AllPredictorCombinations',...
%           'Discriminant','NPredToSample',NPredToSample)
% 
%   T=templateEnsemble(METHOD,NLEARN,LEARNERS,'PARAM1',val1,'PARAM2',val2,...)
%   specifies optional parameter name/value pairs:
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
%                                 pass a scalar, this sampling proportion
%                                 is used for all classes. For example, you
%                                 have class A with 100 observations and
%                                 class B with 10 observations. If you pass
%                                 [2 1] for 'RatioToSmallest', every
%                                 learner in the ensemble is trained on 20
%                                 observations of class A and 10
%                                 observations of class B. If you pass 2 or
%                                 [2 2], every learner is trained on 20
%                                 observations of class A and 20
%                                 observations of class B. Default:
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
%   See also fitcecoc, fitensemble.

%   Copyright 2014 The MathWorks, Inc.

temp = classreg.learning.FitTemplate.make(method,'nlearn',nlearn,'learners',learners,varargin{:});
end
