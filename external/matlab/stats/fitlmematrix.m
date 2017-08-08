function lme = fitlmematrix(X,Y,Z,G,varargin)
%FITLMEMATRIX Create a linear mixed effects model by fitting to data.
%   FITLMEMATRIX is a function for fitting a linear mixed effects model
%   when it is necessary to control the detailed specification of the
%   model. Most models can be fitted more easily using the FITLME function
%   and storing the data in a table array.
%
%   LME = FITLMEMATRIX(X,Y,Z,[]) fits a linear mixed effects model to the N
%   by 1 response vector Y using the N by P fixed effects design matrix X,
%   the N by Q random effects design matrix Z, and returns the fitted model
%   LME. Conceptually, the model is assumed to be:
%
%           Y = X * BETA + Z * B + E   where
%
%        BETA = P by 1 fixed effects vector
%           B = Q by 1 random effects vector 
%           E = N by 1 vector representing the observation error
%
%   Assumptions:
%           (1) The vector B is independent of vector E. 
%           (2) B ~ N(0, SIGMA2 * D)   and 
%           (3) E ~ N(0, SIGMA2 * I)
%   where
%               D = a Q by Q symmetric and positive semi-definite matrix 
%               I = the N by N identity matrix
%          SIGMA2 = the variance of elements of E
%
%   LME = FITLMEMATRIX(X,Y,Z,[]) with Z as a cell array of length R, fits a
%   linear mixed effects model to the N by 1 response vector Y using the N
%   by P fixed effects design matrix X, the N by Q(i) random effects design
%   matrices Z{i}, and returns the fitted model LME. Conceptually, the
%   model is assumed to be:
%
%           Y = X * BETA + SUM_{i = 1 to R} Z{i} * B[i] + E    where
%
%        BETA = P by 1 fixed effects vector 
%        B[i] = a symbol for the Q(i) by 1 random effects vector
%               corresponding to random effects design matrix Z{i}
%           E = N by 1 vector representing the observation error
%
%   Assumptions:
%           (1) The vectors represented by B[i] for i = 1 to R are
%               independent of each other and of the vector E.
%           (2) B[i] ~ N(0, SIGMA2 * D[i])   and 
%           (3)    E ~ N(0, SIGMA2 *  I  )
%   where
%               D[i] = a symbol for a Q(i) by Q(i) symmetric and positive
%                      semi-definite matrix
%                  I = the N by N identity matrix
%             SIGMA2 = the variance of elements of E
%
%   LME = FITLMEMATRIX(X,Y,Z,G) returns a fitted linear mixed effects model
%   LME by fitting to the N by 1 response vector Y using the N by P fixed
%   effects design matrix X and a random effects structure that is
%   specified by both the N by Q random effects design matrix Z and the N
%   by 1 grouping variable G with M "levels" or "groups". G can be a
%   categorical vector, a numeric or logical vector, a char array or a cell
%   array of strings.
%
%   Suppose the symbol L[k] represents the "level" or "group" of grouping
%   variable G for observation k. Conceptually, the response for
%   observation k is modeled as:
%
%           Y(k) = X(k,:) * BETA + Z(k,:) * B[L[k]] + E(k)   where
%
%           BETA = P by 1 fixed effects vector
%        B[L[k]] = a symbol for the Q by 1 random effects vector for level
%                  L[k] of grouping variable G.
%           E(k) = 1 by 1 error term for observation k
%
%   Assumptions:
%           (1) The vectors represented by B[L[k]] for L[k] = level 1 to
%               level M are independent of each other and of the vector E.
%           (2) B[L[k]] ~ N(0, SIGMA2 * D)   and 
%           (3)       E ~ N(0, SIGMA2 * I)
%   where
%               D = a Q by Q symmetric and positive semi-definite matrix 
%               I = the N by N identity matrix
%          SIGMA2 = the variance of elements of E
%
%   LME = FITLMEMATRIX(X,Y,Z,G) with both Z and G as cell arrays of length
%   R, returns a fitted linear mixed effects model LME by fitting to the N
%   by 1 response vector Y using the N by P fixed effects design matrix X.
%   Z is a cell array of length R where Z{j} is a N by Q(j) random effects
%   design matrix. G is a cell array of length R where G{j} is a N by 1
%   grouping variable with M(j) "levels" or "groups". G{j} can be a
%   categorical vector, a numeric vector or a cell array of strings.
%
%   Suppose the symbol L[j,k] represents the "level" or "group" of grouping
%   variable G{j} for observation k. In other words, L[j,k] = G{j}(k).
%   Conceptually, the response vector for the observation k is modeled as:
%
%    Y(k) = X(k,:) * BETA + sum{j = 1 to R} Z{j}(k,:) * B[j, L[j,k]] + E(k)
%   
%    where
%
%            BETA = P by 1 fixed effects vector
%    B[j, L[j,k]] = a Q(j) by 1 random effects vector for grouping variable
%                   G{j} and level L[j,k]. Each level of every grouping
%                   variable gets an independent random effects vector
%            E(k) = 1 by 1 error term for observation k
%
%   Assumptions:
%           (1) The set of random effects vectors B[j, L[j,k]] for j = 1 to
%               R and L[j,k] = level 1 to level M(j) are independent of
%               each other and of the vector E.
%           (2) B[j, L[j,k]] ~ N( 0, SIGMA2 * D[j] )   and 
%           (3)            E ~ N( 0, SIGMA2 *   I  )
%   where
%           D[j] = a symbol for the Q(j) by Q(j) symmetric and positive
%                  semi-definite matrix
%              I = N by N identity matrix
%         SIGMA2 = the variance of elements of E
%
%   LME = FITLMEMATRIX(X,Y,Z,G, PARAM1,VAL1, PARAM2,VAL2,...) specifies one
%   or more of the following name/value pairs:
%
%       'FixedEffectPredictors'    
%                         A cell array of length P containing the names of
%                         the columns of the fixed effects design matrix X.
%                         Default is {'x1','x2',...,'xP'}.
%
%       'RandomEffectPredictors'    
%                         If Z is a single N by Q design matrix,
%                         'RandomEffectPredictors' is a cell array of
%                         length Q containing the names of the columns of
%                         Z. Default is {'z1',...,'zQ'}. If Z is a cell
%                         array of length R with Z{i} a N by Q(i) design
%                         matrix, then 'RandomEffectPredictors' is a cell
%                         array of length R, element i of which is a
%                         cell array containing the names of the columns of
%                         Z{i}. Default is illustrated by example: If Q(1)
%                         = 2 and Q(2) = 3 then the default names for
%                         columns of Z{1} are: {'z11','z12'} and for
%                         columns of Z{2} are: {'z21','z22','z23'}. 
%   
%       'ResponseVarName'     
%                         A name for the response variable Y. Default is 'y'. 
%
%       'RandomEffectGroups'  
%                         When G is a vector, 'RandomEffectGroups' is a 
%                         single string naming G. Default is 'g'. When G is
%                         a cell array of length R, 'RandomEffectGroups' is
%                         a cell array of length R, element i of which is a
%                         name for the grouping variable G{i}. Default is
%                         {'g1',...,'gR'}.
%
%       'CovariancePattern'   
%                         If Z is a cell array of length R with Z{i} a N by
%                         Q(i) design matrix, 'CovariancePattern' is a cell
%                         array of length R such that element i of this
%                         cell array specifies the pattern of the
%                         covariance matrix D[i] of the random effects
%                         vector associated with Z{i}. Each element of the
%                         cell array for 'CovariancePattern' can either be
%                         a string or a logical matrix. Allowed values for
%                         element i are:
%
%           'FullCholesky'  - a full covariance matrix using the Cholesky 
%                             parameterization (Default). All elements of 
%                             the covariance matrix are estimated.
%
%           'Full'          - a full covariance matrix using the 
%                             log-Cholesky parameterization. All elements 
%                             of the covariance matrix are estimated.
%
%           'Diagonal'      - a diagonal covariance matrix. Off-diagonal 
%                             elements of the covariance matrix are 
%                             constrained to be 0.
%
%           'Isotropic'     - a diagonal covariance with equal variances.
%                             Off-diagonal elements are constrained to be 
%                             0 and diagonal elements are constrained to be 
%                             equal.
%
%           'CompSymm'      - a compound symmetry structure i.e., common 
%                             variance along diagonals and equal 
%                             correlation between all random effects.
%
%           PAT             - A square symmetric logical matrix. If 
%                             PAT(a,b) = false then the (a,b) element of 
%                             the covariance matrix is constrained to 0.
%
%                         When Z is a N by Q matrix, 'CovariancePattern'
%                         can be a string or a logical matrix instead of a
%                         1 by 1 cell array.
%
%       'FitMethod'       Specifies the method to use for estimating linear
%                         mixed effects model parameters. Choices are:
%
%           'ML'   - maximum likelihood (Default) 
%
%           'REML' - restricted maximum likelihood
%   
%       'Weights'         Vector of N non-negative weights, where N is the
%                         number of rows in Y. Default is ones(N,1).
%
%       'Exclude'         Vector of integer or logical indices into the
%                         rows of X, Y and Z, that should be excluded from
%                         the fit. Default is to use all rows.
%
%       'DummyVarCoding'  A string specifying the coding to use for dummy
%                         variables created from categorical variables.
%                         Valid coding schemes are 'reference' (coefficient
%                         for first category set to zero), 'effects'
%                         (coefficients sum to zero) and 'full' (one dummy
%                         variable for each category). Default is
%                         'reference'.
%
%       'Optimizer'       A string specifying the algorithm to use for
%                         optimization. Valid values of this parameter are
%                         'quasinewton' (Default) and 'fminunc'. If
%                         'Optimizer' is 'quasinewton', a trust region
%                         based quasi-Newton optimizer is used. If you have
%                         the Optimization Toolbox, you can also specify
%                         that fminunc be used for optimization by setting
%                         'Optimizer' value to 'fminunc'.
%
%       'OptimizerOptions'
%                         If 'Optimizer' is 'quasinewton', then
%                         'OptimizerOptions' is a structure created by
%                         statset('fitlmematrix'). The quasi-Newton
%                         optimizer uses the following fields:
%
%           'TolFun'        - Relative tolerance on the gradient of the 
%                             objective function. Default is 1e-6.
%                  
%           'TolX'          - Absolute tolerance on the step size.
%                             Default is 1e-12.
%
%           'MaxIter'       - Maximum number of iterations allowed. 
%                             Default is 10000.
%
%           'Display'       - Level of display.  'off', 'iter', or 'final'.
%                             Default is off.
%
%                         If 'Optimizer' is 'fminunc', then 
%                         'OptimizerOptions' is an object set up using 
%                         optimoptions('fminunc'). See the documentation
%                         for optimoptions for a list of all the options
%                         supported by fminunc.
%
%                         If 'OptimizerOptions' is not supplied and
%                         'Optimizer' is 'quasinewton' then the default
%                         options created by statset('fitlmematrix') are
%                         used. If 'OptimizerOptions' is not supplied and
%                         'Optimizer' is 'fminunc' then the default options
%                         created by optimoptions('fminunc') are used with
%                         the 'Algorithm' set to 'quasi-newton'.
%
%       'StartMethod'     Method to use to start iterative optimization.
%                         Choices are 'default' (Default) and 'random'. If
%                         'StartMethod' is 'random', a random initial value
%                         is used to start iterative optimization,
%                         otherwise an internally defined default value is
%                         used.
%
%       'Verbose'         Either true or false. If true, then iterative
%                         progress of the optimizer is displayed on screen.
%                         The setting for 'Verbose' overrides field
%                         'Display' in 'OptimizerOptions'. Default is
%                         false.
%
%      'CheckHessian'     Either true or false. If 'CheckHessian' is true
%                         then optimality of the solution is verified by
%                         performing positive definiteness checks on the
%                         Hessian of the objective function with respect to
%                         unconstrained parameters at convergence. Hessian
%                         checks are also performed to determine if the
%                         model is overparameterized in the number of
%                         covariance parameters. Default is false.
%
%   Example: Three ways to specify equivalent models
%       % Fit model from dataset using fitlme
%       load carsmall
%       ds = dataset(MPG, Model_Year, Weight);
%       lme1 = fitlme(ds,'MPG~Weight+(1|Model_Year)')
% 
%       % Fit the same model by specifying data as matrices
%       y = MPG;
%       x = [ones(size(Weight)), Weight];
%       z = ones(size(y));
%       lme2 = fitlmematrix(x,y,z,Model_Year)
% 
%       % Fit the same model by building the grouping into the Z matrix
%       z = double([Model_Year==70, Model_Year==76, Model_Year==82]);
%       lme3 = fitlmematrix(x,y,z,[],'covariancepattern','isotropic')
%
%   See also FITLME, LinearMixedModel.

%   Copyright 2012-2014 The MathWorks, Inc.

    narginchk(4,Inf);
    lme = LinearMixedModel.fitmatrix(X,Y,Z,G,varargin{:});

end % end of fitlmematrix.
