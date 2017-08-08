function model = stepwiselm(X,varargin) % [X, y | DS], start, ...
%STEPWISELM Create linear regression model by stepwise regression.
%   LM = STEPWISELM(DS,MODELSPEC) fits the model specified by MODELSPEC to
%   variables in the dataset/table DS, adds or removes terms by stepwise
%   regression, and returns the linear model LM. MODELSPEC can be any of
%   the values accepted by the FITLM function. The default is 'constant' to
%   start with only the constant term in the model.
%
%   After the initial fit, the function examines a set of available terms,
%   and adds the best one to the model if an F-test for adding the term has
%   a p-value 0.05 or less. If no terms can be added, it examines the terms
%   currently in the model, and removes the worst one if an F-test for
%   removing it has a p-value 0.10 or greater. It repeats this process
%   until no terms can be added or removed. The function never removes the
%   constant term. It may add terms defined by MODELSPEC, as well as terms
%   that are two-way interactions among the predictors.
%
%   LM = STEPWISELM(X,Y) fits a linear regression model using the column
%   vector Y as a response variable and the columns of the matrix X as
%   predictor variables, performs stepwise regression, and returns the
%   final result as the linear model LM.
%
%   LM = STEPWISELM(X,Y,MODELSPEC) uses the model specified by MODELSPEC as
%   the initial model. See FITLM for valid MODELSPEC values.
%
%   LM = STEPWISELM(...,PARAM1,VAL1,PARAM2,VAL2,...) specifies one or more
%   of the following name/value pairs:
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
%                         or remove, chosen from 'SSE' (default), 'AIC',
%                         'BIC', 'Rsquared', 'AdjRsquared'.
%      'PEnter'           For the 'SSE' criterion, a value E such that a
%                         term may be added if its p-value is less than or
%                         equal to E. For the other criteria, a term may be
%                         added if the improvement in the criterion is at
%                         least E.
%      'PRemove'          For the 'SSE' criterion, a value R such that a
%                         term may be removed if its p-value is greater or
%                         equal to R. For the other criteria, a term may be
%                         added if it reduces the criterion no more than R.
%      'NSteps'           The maximum number of steps that may be taken,
%                         starting from the initial model. Default is no
%                         limit on the number of steps.
%      'Verbose'          An integer from 0 to 2 controlling the display of
%                         information. Verbose=1 (the default) displays the
%                         action taken at each step. Verbose=2 also
%                         displays the actions evaluated at each step.
%                         Verbose=0 suppresses all display.
%
%   The following table shows the default 'PEnter' and 'PRemove' values for
%   the different criteria, and indicates which must be larger than the
%   other:
%
%      Criterion     PEnter   PRemove    Compared against
%      'SSE'         0.05   < 0.10       p-value for F test
%      'AIC'         0      < 0.01       change in AIC
%      'BIC'         0      < 0.01       change in BIC
%      'Rsquared'    0.1    > 0.05       increase in R-squared
%      'AdjRsquared' 0      > -0.05      increase in adjusted R-squared
%
%   Example:
%      % Perform stepwise regression, starting from the constant model
%      % (intercept only), adding linear terms as required.
%      load hald
%      lm = stepwiselm(ingredients,heat,'constant','upper','linear')
%
%   See also fitlm, LinearModel.

%   Copyright 2011-2014 The MathWorks, Inc.

model = LinearModel.stepwise(X,varargin{:});
end