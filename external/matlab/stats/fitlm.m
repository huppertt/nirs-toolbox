function model = fitlm(X,varargin)
%FITLM Create linear regression model by fitting to data.
%   LM = FITLM(DS,MODELSPEC) fits the model specified by MODELSPEC to
%   variables in the dataset/table array DS, and returns the linear model LM.
%   MODELSPEC can be any of the following:
%
%       'linear'            Linear (main effect) terms only.
%       'interactions'      Linear and pairwise interaction terms.
%       'purequadratic'     Linear and squared terms.
%       'quadratic'         Linear, pairwise interactions, and squares.
%       'polyIJK...'        Polynomial with all terms up to power I for the
%                           first predictor, J for the second, K for the
%                           third, and so on.
%       FORMULA             a string such as 'y ~ x1 + x2 + x3*x4' defining
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
%       TERMS               A T-by-V matrix where T is the desired number
%                           of terms in the model and V is the number of
%                           variables in DS. The (I,J) element indicates
%                           the power of variable J in term I.
%
%   The following are some examples of FORMULA term expressions when the
%   predictors are x1, x2, and x3:
%       'x1+x2+x3'             Linear model including constant
%       'x1+x2+x3-1'           Linear model without constant
%       'x1^2+x2^2'            Constant, linear, squared terms
%       'x1*x2*x3-x1:x2:x3'    All except the three-way interaction
%       'x1*(x2+x3)'           Linear, plus two-way interactions with x1
%
%   The following are some examples of TERMS matrix rows when the dataset/table
%   DS consists of predictors x1, x2, and x3, followed by response y:
%       [0 0 0 0] represents a constant term or intercept
%       [0 1 0 0] represents x2
%       [1 0 1 0] represents the product x1:x3
%       [2 0 0 0] represents x1:x1
%       [0 1 2 0] represents x2:x3:x3
%
%   LM = FITLM(X,Y) fits a linear regression model using the column vector
%   Y as a response variable and the columns of the matrix X as predictor
%   variables, and returns the linear model LM.
%
%   LM = FITLM(X,Y,MODELSPEC) fits the model specified by MODELSPEC. If
%   MODELSPEC is a formula, it should specify the response as y and the
%   predictors as x1,x2,.... The default is 'linear'.
%
%   LM = FITLM(...,PARAM1,VAL1,PARAM2,VAL2,...) specifies one or more of
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
%      'RobustOpts'       'off' (default) to perform least squares, 'on' to
%                         perform robust regression using the bisquare
%                         weighting function, or a structure.
%
%   When the 'RobustOpts' value is a structure, it must have a
%   'RobustWgtFun' field and may have a 'Tune' field. The value of the
%   'RobustWgtFun' field can be any of 'andrews', 'bisquare', 'cauchy',
%   'fair', 'huber', 'logistic', 'talwar', or 'welsch', or it can be a
%   function that takes a residual vector as input and produces a weight
%   vector as output.  The value of the 'Tune' field is a tuning constant.
%   The residuals are scaled by the tuning constant (default depends on
%   weighting function and is 1 for a function handle) and by an
%   estimate of the error standard deviation before the weight function is
%   called.  'RobustWgtFun' can be specified using @ (as in @myfun).
%
%   Examples:
%      % Fit to data in matrices
%      load hald
%      lm = fitlm(ingredients,heat)
%
%      % Fit to data in a dataset
%      load carsmall
%      d = dataset(MPG,Weight);
%      d.Year = ordinal(Model_Year);
%      lm = fitlm(d,'MPG ~ Year + Weight + Weight^2')
%
%      % Fit to data in a table
%      load carsmall
%      d = table(MPG,Weight);
%      d.Year = ordinal(Model_Year);
%      lm = fitlm(d,'MPG ~ Year + Weight + Weight^2')
%
%   See also LinearModel, STEPWISELM.

%   Copyright 2013-2014 The MathWorks, Inc.

model = LinearModel.fit(X,varargin{:});
end