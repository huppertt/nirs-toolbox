function model = stepwiseglm(X,varargin) % [X, y | DS], start, ...
%STEPWISEGLM Create generalized linear regression model by stepwise regression.
%   GLM = STEPWISEGLM(DS,MODELSPEC) fits the model specified by MODELSPEC
%   to variables in the dataset/table DS, adds or removes terms by stepwise
%   regression, and returns the generalized linear model GLM. MODELSPEC can
%   be any of the values accepted by the FITGLM function. The default is
%   'constant' to start with only the constant term in the model.
%
%   After the initial fit, the function examines a set of available terms,
%   and adds the best one to the model if a test for adding the term has
%   a p-value 0.05 or less. If no terms can be added, it examines the terms
%   currently in the model, and removes the worst one if a test for
%   removing it has a p-value 0.10 or greater. It repeats this process
%   until no terms can be added or removed. The function never removes the
%   constant term. It may add terms defined by MODELSPEC, as well as terms
%   that are two-way interactions among the predictors.
%
%   GLM = STEPWISEGLM(X,Y) fits a linear regression model using the column
%   vector Y as a response variable and the columns of the matrix X as
%   predictor variables, performs stepwise regression, and returns the
%   final result as the generalized linear model GLM.
%
%   GLM = STEPWISEGLM(X,Y,MODELSPEC) uses the model specified by MODELSPEC
%   as the initial model. See FITGLM for valid MODELSPEC values.
%
%   GLM = STEPWISEGLM(...,MODELSPEC,PARAM1,VAL1,PARAM2,VAL2,...) specifies
%   one or more of the following name/value pairs:
%
%      'Weights'          Vector of N non-negative weights, where N is the
%                         number of rows in DS or Y. Default is ones(N,1).
%      'VarNames'         Cell array of strings specifying the names to use
%                         for the columns in X. Default is {'x1','x2',...}
%                         for the predictors and 'y' for the response.
%                         Not allowed when fitting to a dataset/table.
%      'CategoricalVars'  Vector of integer or logical indices specifying
%                         the variables in DS or the columns in X that
%                         should be treated as categorical. Default is to
%                         treat DS variables as categorical if they are
%                         categorical, logical, or char arrays, or cell
%                         arrays of strings.
%      'Exclude'          Vector of integer or logical indices into the
%                         rows of DS, or X and Y, that should be excluded
%                         from the fit. Default is to use all rows.
%      'Intercept'        true (default) to include a constant term in the
%                         model, or false to omit it.
%      'PredictorVars'    A specification of the variables to use as
%                         predictors, either as a cell array of variable
%                         names, or a vector of integer or logical indices
%                         into the variables in DS or the columns in X.
%      'ResponseVar'      The response variable, specified either as a
%                         variable name or number.
%      'Lower'            A model specification defining the terms that
%                         cannot be removed from the model. Default
%                         'constant', meaning only the intercept.
%      'Upper'            A model specification defining the terms
%                         available to be added to the model. Default
%                         'interactions' for pairwise interaction terms.
%      'Criterion'        The criterion to be used in choosing terms to add
%                         or remove, chosen from 'Deviance' (default),
%                         'SSE', 'AIC', 'BIC', 'Rsquared', 'AdjRsquared'.
%      'PEnter'           For the 'Deviance' or 'SSE' criterion, a value E
%                         such that a term may be added if its p-value is
%                         less than or equal to E. For the other criteria,
%                         a term may be added if the improvement in the
%                         criterion is at least E.
%      'PRemove'          For the 'Deviance' or 'SSE' criterion, a value R
%                         such that a term may be removed if its p-value is
%                         greater or equal to R. For the other criteria, a
%                         term may be added if it reduces the criterion no
%                         more than R.
%      'NSteps'           The maximum number of steps that may be taken,
%                         starting from the initial model. Default is no
%                         limit on the number of steps.
%      'Verbose'          An integer from 0 to 2 controlling the display of
%                         information. Verbose=1 (the default) displays the
%                         action taken at each step. Verbose=2 also
%                         displays the actions evaluated at each step.
%                         Verbose=0 suppresses all display.
%      'Distribution'     Name of the distribution of the response, chosen
%                         from the following:
%                 'normal'             Normal distribution (default)
%                 'binomial'           Binomial distribution
%                 'poisson'            Poisson distribution
%                 'gamma'              Gamma distribution
%                 'inverse gaussian'   Inverse Gaussian distribution
%      'BinomialSize'     Vector or name of a variable of the same length
%                         as the response, specifying the size of the sample
%                         (number of trials) used in computing Y. This is the
%                         N parameter for the fitted binomial distribution
%                         and is accepted only when the 'Distribution'
%                         parameter is 'binomial'. May be a scalar if all
%                         observations have the same value. Default is 1. As
%                         an alternative to this parameter, you can specify
%                         the response as a two-column vector with counts
%                         in column 1 and BinomialSize in column 2.
%      'Link'             The link function to use in place of the
%                         canonical link. The link function defines the
%                         relationship f(mu) = x*b between the mean
%                         response mu and the linear combination of
%                         predictors x*b.  Specify the link as:
%                 'identity'    f(mu) = mu
%                 'log'         f(mu) = log(mu)
%                 'logit'       f(mu) = log(mu/(1-mu))
%                 'probit'      f(mu) = norminv(mu)
%                 'comploglog'  f(mu) = log(-log(1-mu))
%                 'reciprocal'  f(mu) = log(-log(mu))
%                 P (number)    f(mu) = mu.^P
%                 S (struct)    structure with three fields whose values
%                               are function handles and with these names:
%                                  S.Link           link function
%                                  S.Derivative     derivative
%                                  S.Inverse        inverse of link
%      'Offset'           Vector or name of a variable with the same
%                         length as the response. This is used as an
%                         additional predictor with a coefficient value
%                         fixed at 1.0.
%      'DispersionFlag'   true to estimate a dispersion parameter when
%                         computing standard errors in models with the
%                         binomial and Poisson distributions, or false
%                         (default) to use the theoretical value. The
%                         function always estimates the dispersion for
%                         other distributions.
%
%   The following table shows the default 'PEnter' and 'PRemove' values for
%   the different criteria, and indicates which must be larger than the
%   other:
%
%      Criterion     PEnter   PRemove    Compared against
%      'Deviance'    0.05   < 0.10       p-value for F or chi-square test
%      'SSE'         0.05   < 0.10       p-value for F test
%      'AIC'         0      < 0.01       change in AIC
%      'BIC'         0      < 0.01       change in BIC
%      'Rsquared'    0.1    > 0.05       increase in R-squared
%      'AdjRsquared' 0      > -0.05      increase in adjusted R-squared
%
%   Example:
%      % See if we can identify the correct predictors using random data
%      % with a Poisson regression model.
%      x = randn(100,20);
%      mu = exp(x(:,[5 10 15])*[.4;.2;.3] + 1);
%      y = poissrnd(mu);
%      stepwiseglm(x,y,'constant','upper','linear','distr','poisson')
%
%   See also FITGLM, GeneralizedLinearModel.

%   Copyright 2013-2014 The MathWorks, Inc.

model = GeneralizedLinearModel.stepwise(X,varargin{:});
end