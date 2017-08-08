function model = fitglm(X,varargin) % [X, y | DS], modelDef, ...
%FITGLM Create generalized linear regression model by fitting to data.
%   GLM = FITGLM(DS,MODELSPEC) fits the model specified by MODELSPEC to
%   variables in the dataset/table array DS, and returns the generalized linear
%   model GLM. The fit uses a normal distribution and an identity link
%   function. MODELSPEC can be any of the following:
%
%       'linear'            Linear (main effect) terms only.
%       'interactions'      Linear and pairwise interaction terms.
%       'purequadratic'     Linear and squared terms.
%       'quadratic'         Linear, pairwise interactions, and squares.
%       'polyIJK...'        Polynomial with all terms up to power I for the
%                           first predictor, J for the second, K for the
%                           third, and so on.
%       formula             a string such as 'y ~ x1 + x2 + x3*x4' defining
%                           the response and the predictor terms. A formula
%                           string always has the response variable name,
%                           followed by '~', followed by one or more terms
%                           joined by '+' or '-'. The following are the
%                           rules for constructing a formula:
%                             A + B       term A and term B
%                             A - B       term A but without term B
%                             A:B         the product of A and B
%                             A*B         A + B + A:B
%                             A^2         A + A:A
%                             ()          grouping of terms
%
%   The following are some examples of formulas when the predictors are x1,
%   x2, and x3:
%       'x1+x2+x3'             Linear model including constant
%       'x1+x2+x3-1'           Linear model without constant
%       'x1^2+x2^2'            Constant, linear, squared terms
%       'x1*x2*x3-x1:x2:x3'    All except the three-way interaction
%       'x1*(x2+x3)'           Linear, plus two-way interactions with x1
%
%   GLM = FITGLM(X,Y) fits a generalized linear regression model using the
%   column vector Y as a response variable and the columns of the matrix X
%   as predictor variables, and returns the generalized linear model GLM.
%
%   GLM = FITGLM(X,Y,MODELSPEC) fits the model specified by MODELSPEC. If
%   MODELSPEC is a formula, it should specify the response as y and the
%   predictors as x1,x2,.... The default is 'linear'.
%
%   GLM =FITGLM(...,PARAM1,VAL1,PARAM2,VAL2,...) specifies one or more of
%   the following name/value pairs:
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
%                 'loglog'      f(mu) = log(-log(mu))
%                 'reciprocal'  f(mu) = mu.^(-1)
%                 P (number)    f(mu) = mu.^P
%                 S (struct)    structure with three fields whose values
%                               are function handles and with these names:
%                                  S.Link           link function
%                                  S.Derivative     derivative
%                                  S.Inverse        inverse of link
%                         The default is the canonical link that depends on
%                         the distribution:
%                 'identity'    normal distribution
%                 'logit'       binomial distribution
%                 'log'         Poisson distribution
%                 -1            gamma distribution
%                 -2            inverse gaussian distribution
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
%   Example:
%      % Fit the probability of smoking using a logistic model with
%      % all except the three-way interaction of three predictors.
%      load hospital
%      formula = 'Smoker ~ Age*Weight*Sex - Age:Weight:Sex';
%      glm = fitglm(hospital,formula,'distr','binomial')
%
%   See also GeneralizedLinearModel, STEPWISEGLM.

%   Copyright 2013-2014 The MathWorks, Inc.

model = GeneralizedLinearModel.fit(X,varargin{:});
end