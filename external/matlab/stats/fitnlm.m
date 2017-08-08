function model = fitnlm(X,varargin) % [X, y | DS], modelDef, initialCoefs, ...
%FITNLM Create nonlinear regression model by fitting to data.
%   NLM = FITNLM(DS,MODELFUN,BETA0) fits the model specified by MODELFUN to
%   variables in the dataset/table array DS, and returns the nonlinear model NLM.
%   The coefficients are estimated using an iterative procedure starting
%   from the initial values in the vector BETA0.
%
%   MODELFUN can be either of the following:
%
%     1. A function, specified using @, that accepts two arguments, a
%        coefficient vector and an array of predictor values, and computes
%        a vector of Y values. 
%     2. A text string, such as 'y ~ b0+b1*sin(b2*X)', defining the
%        response and a mathematical expression for the function. Symbols
%        in the function expression that match variable names in the
%        dataset/table are taken to be predictors, and the others are
%        coefficients to be estimated.
%
%   NLM = FITNLM(X,Y,MODELFUN,BETA0) fits a nonlinear regression model
%   using the column vector Y as a response variable and the columns of the
%   matrix X as predictor variables.
%
%   NLM = FITNLM(...,BETA0,PARAM1,VAL1,PARAM2,VAL2,...) specifies one or
%   more of the following name/value pairs:
%
%      'CoefficientNames' Cell array of strings specifying the names of
%                         the coefficients.
%      'Weights'          Vector W of N non-negative weights, where N is the
%                         number of rows in DS or Y. Default is ones(N,1).
%                         The error variance at observation i is estimated as 
%                         MSE * (1/W(i)) where MSE is the mean squared error.  
%                         Weights can also be specified as a function handle 
%                         that accepts a vector of predicted response values 
%                         and returns a vector of real positive weights as 
%                         output. 
%      'VarNames'         Cell array of strings specifying the names to use
%                         for the columns in X. Default is {'x1','x2',...}
%                         for the predictors and 'y' for the response.
%                         Not allowed when fitting to a dataset/table.
%      'Exclude'          Vector of integer or logical indices into the
%                         rows of DS, or X and Y, that should be excluded
%                         from the fit. Default is to use all rows.
%      'PredictorVars'    A specification of the variables to use as
%                         predictors, either as a cell array of variable
%                         names, or a vector of integer or logical indices
%                         into the variables in DS or the columns in X.
%      'ResponseVar'      The response variable, specified either as a
%                         variable name or number.
%      'Options'          A structure created by STATSET to control the
%                         iterative fitting procedure including iteration
%                         limits, convergence tolerance, and robust
%                         fitting. The fields relevant here are the ones
%                         created by calling statset('fitnlm').
%      'ErrorModel'       A string specifying the form of the error term.
%                         Default is 'constant'. Each model defines the error
%                         using a standard zero mean and unit variance variable e
%                         with independent components, function value f, and one 
%                         or two parameters a and b. The Choices are:              
%                    
%                         'constant' (default)         y = f + a*e
%                         'proportional'               y = f + b*f*e
%                         'combined'                   y = f + (a+b*abs(f))*e
%
%                         If not specified, a 'constant' error model will be 
%                         used. The only allowed 'ErrorModel' when using 'Weights' 
%                         is 'constant'. OPTIONS.RobustWgtFun must be [] when 
%                         using error models other than 'constant'.
%   'ErrorParameters'     A numeric array containing initial estimates of the
%                         error model parameters of the chosen 'ErrorModel'. 
%                         Specify this as:   
%                     
%                         a       (default 1)       for 'constant'
%                         b       (default 1)       for 'proportional'
%                         [a b]   (default [1,1])   for 'combined'
% 
%                         For example, if 'ErrorModel' is 'combined', one could 
%                         use 'ErrorParameters' as [1,2]. If not specified, the 
%                         default values mentioned in brackets above will be used 
%                         as initial estimates.
%
%   Example:
%      % Fit the Hougen model to reaction data
%      load reaction
%      nlm = fitnlm(reactants,rate,@hougen,[1 .05 .02 .1 2])
%
%      % Same fit using a text representation of the model
%      myfun = 'rate~(b1*x2-x3/b5)/(1+b2*x1+b3*x2+b4*x3)';
%      nlm = fitnlm(reactants,rate,myfun,[1 .05 .02 .1 2])
%
%   See also NonLinearModel, statset.

%   Copyright 2013-2014 The MathWorks, Inc.

model = NonLinearModel.fit(X,varargin{:});
end
