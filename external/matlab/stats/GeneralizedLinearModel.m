classdef (Sealed = true) GeneralizedLinearModel < classreg.regr.TermsRegression
%GeneralizedLinearModel Fitted multiple generalized linear regression model.
%   GLM = FITGLM(...) fits a generalized linear model to data. The fitted
%   model GLM is a GeneralizedLinearModel that can predict a response as a
%   function of predictor variables and terms created from predictor
%   variables. Generalized linear models support several distributions for
%   the response, and several choices of link function for relating the
%   predictors to the parameters of the response distribution.
%
%   GeneralizedLinearModel methods:
%       addTerms - Add terms to linear model
%       removeTerms - Remove terms from linear model
%       step - Selectively add or remove terms from linear model
%       coefCI - Coefficient confidence intervals
%       coefTest - Linear hypothesis test on coefficients
%       predict - Compute predicted values given predictor values
%       feval - Evaluate model as a function
%       random - Generate random response values given predictor values
%       plotDiagnostics - Plot diagnostics of fitted model
%       plotResiduals - Plot residuals of fitted model
%       plotSlice - Plot slices through the fitted regression surface
%
%   GeneralizedLinearModel properties:
%       Coefficients - Coefficients and related statistics
%       Rsquared - R-squared and adjusted R-squared
%       ModelCriterion - AIC and other model criteria
%       Fitted - Vector of fitted (predicted) values
%       Residuals - Table containing residuals of various types
%       ResponseName - Name of response
%       PredictorNames - Names of predictors
%       NumPredictors - Number of predictors
%       Variables - Table of variables used in fit
%       NumVariables - Number of variables used in fit
%       VariableNames - Names of variables used in fit
%       VariableInfo - Information about variables used in the fit
%       NumObservations - Number of observations in the fit
%       ObservationNames - Names of observations in the fit
%       ObservationInfo - Information about observations used in the fit
%       DFE - Degrees of freedom for error
%       SSE - Error sum of squares
%       SST - Total sum of squares
%       SSR - Regression sum of squares
%       Steps - Stepwise fit results
%       Formula - Representation of the model used in this fit
%       LogLikelihood - Log of likelihood function at coefficient estimates
%       CoefficientCovariance - Covariance matrix for coefficient estimates
%       CoefficientNames - Coefficient names
%       NumCoefficients - Number of coefficients
%       NumEstimatedCoefficients - Number of estimated coefficients
%       Distribution - Distribution of the response
%       Link - Link relating the distribution parameters to the predictors
%       Dispersion - Theoretical or estimated dispersion parameter
%       DispersionEstimated - Flag indicating if Dispersion was estimated
%       Deviance - Deviance of the fit
%
%   See also FITGLM, STEPWISEGLM, LinearModel, NonLinearModel.

%   Copyright 2011-2013 The MathWorks, Inc.

    properties(Constant,Hidden)
        SupportedResidualTypes = {'Raw' 'LinearPredictor' 'Pearson' 'Anscombe' 'Deviance'};
    end

    properties(GetAccess='public',SetAccess='protected')
%Dispersion Parameter defining the variance of the response.
%    The value of the Dispersion property depends on the Distribution of
%    the model. For a normal distribution, the Dispersion is the mean
%    squared error of the residuals. For other distributions, the
%    Dispersion multiplies the variance function for the distribution.
%
%    For example, the variance function for the binomial distribution is
%    P*(1-P)/N where P is the probability parameter and N is the sample
%    size parameter. If the Dispersion parameter is near 1, then the
%    variance of the data appears to agree with the theoretical variance of
%    the binomial distribution. If the Dispersion parameter is larger than
%    1, the data are said to be overdispersed relative to the binomial
%    distribution.
%
%    See also GeneralizedLinearModel, DispersionEstimated.
        Dispersion = 0;

%DispersionEstimated Estimated dispersion used to estimate standard errors.
%    The DispersionEstimated property is true if the Dispersion parameter
%    is used in computing standard errors for the coefficients, or false if
%    it is not. This property can be false only for the binomial and
%    Poisson distributions.
%
%    See also GeneralizedLinearModel, Dispersion.
        DispersionEstimated = false;
        
%Deviance Deviance of the fit.
%    The Deviance property is the deviance of the fit, which is equal to
%    -2 times the log likelihood.
%
%    The deviance is useful for comparing two models when one is a special
%    case of the other. The difference D between the deviance of the two
%    models is -2 times the log of the likelihood ratio. Asymptotically, D
%    has a chi-square distribution with degrees of freedom V equal to the
%    number of parameters that are estimated in one model but fixed
%    (typically at 0) in the other. The p-value for this test is
%    1-chi2cdf(D,V).
%
%    See also GeneralizedLinearModel, NumEstimatedCoefficients, chi2cdf.
        Deviance = NaN;

%Offset Predictor value with coefficient known to be 1.
%    The Offset property is a vector specifying the values of a special
%    predictor whose coefficient is not estimated, but is set to the value
%    1.
%
%    For example, consider a Poisson regression model. Suppose the number
%    of counts is known for theoretical reasons to be proportional to a
%    predictor A. By using the log link function and by specifying log(A)
%    as an offset, you can force the model to satisfy this theoretical
%    constraint.
%
%    See also GeneralizedLinearModel, Link.
        Offset = 0;
    end
    properties(GetAccess='protected',SetAccess='protected')
        DistributionName = 'normal';
        DevianceNull = NaN;
        IRLSWeights = []; % not reduced with Subset
    end
    properties(Dependent,GetAccess='public',SetAccess='protected')
%Distribution Response distribution and related information.
%    The Distribution property is a structure providing the name and other
%    characteristics of the distribution of the response. This structure
%    has three fields:
%
%        Name             Name of the distribution; one of 'normal',
%                         'binomial', 'poisson', 'gamma', or 'inverse gamma'.
%        DevianceFunction Function that computes the components of the
%                         deviance as a function of the fitted parameter
%                         values and the response values.
%        VarianceFunction Function that computes the theoretical variance
%                         for the distribution as a function of the fitted
%                         parameter values. When the DispersionEstimated
%                         parameter is true, the Dispersion parameter
%                         multiplies the variance function in the
%                         computation of the coefficient standard errors.        
%
%    See also GeneralizedLinearModel, Link, Dispersion, DispersionEstimated.
        Distribution

%Link Link between the distribution parameter and predictor values.
%    The Link property is a structure providing the name and other
%    characteristics of the link function. The link is a function F that
%    links the distribution parameter MU to the fitted linear combination
%    XB of the predictors:  F(MU) = XB.  The structure has four fields:
%
%        Name             Name of the link function, or '' if the link is
%                         specified as functions rather than by name.
%        LinkFunction     The function that defines F.
%        DevianceFunction Derivative of F.
%        VarianceFunction Inverse of F.
%
%    See also GeneralizedLinearModel, Distribution.
        Link
        
%Residuals - Residual values.
%   The Residuals property is a table of residuals. It is a table that
%   has one row for each observation and the following variables:
%
%      'Raw'             Observed minus fitted values
%      'LinearPredictor' Residuals on the linear predictor scale, equal to
%                        the adjusted response value minus the fitted linear
%                        combination of the predictors.
%      'Pearson'         Residuals divided by the estimated standard deviation
%                        of the response.
%      'Anscombe'        Residuals defined on transformed data with the
%                        transformation chosen to remove skewness.
%      'Deviance'        Residuals based on the contribution of each
%                        observation to the deviance.
%
%   To obtain any of these columns as a vector, index into the property
%   using dot notation. For example, in the model M, the ordinary or
%   raw residual vector is
%
%      r = M.Residuals.Raw
%
%   See also GeneralizedLinearModel, plotResiduals, Fitted, predict, random.
        Residuals

%Fitted - Fitted (predicted) values.
%   The Fitted property is a table of fitted values. It is a table that
%   has one row for each observation and the following variables:
%
%      'Response'         Fitted values on the scale of the response.
%      'LinearPredictor'  Fitted values on the scale of the linear
%                         predictor. These are the same as the link
%                         function applied to the 'Response' fitted values.
%      'Probability'      Fitted probabilities (this column is included
%                         only with the binomial distribution).
%
%   To obtain any of the columns as a vector, index into the property
%   using dot notation. For example, in the model M, the vector of fitted
%   values on the response scale is
%
%      f = M.Fitted.Response
%
%   The fitted values are computed using the predictor values used to fit
%   the model. Use the PREDICT method to compute predictions for other
%   predictor values and to compute confidence bounds for the predicted
%   values.
%
%   See also GeneralizedLinearModel, Residuals, predict, random.
        Fitted
        
%Diagnostics - Regression diagnostics.
%   The Diagnostics property is a structure containing a set of diagnostics
%   helpful in finding outliers and influential observations. The structure
%   contains the following fields:
%      Leverage  Diagonal elements of the Hat matrix
%      CooksDistance Cook's measure of scaled change in fitted values
%      HatMatrix Projection matrix to compute fitted from observed responses
%
%   Leverage indicates to what extent the predicted value for an
%   observation is determined by the observed value for that observation. A
%   value close to 1 indicates that the prediction is largely determined by
%   that observation, with little contribution from the other observations.
%   A value close to 0 indicates the fit is largely determined by the other
%   observations. For a model with P coefficients and N observations, the
%   average value of Leverage is P/N. Observations with Leverage larger
%   than 2*P/N may be considered to have high leverage.
%
%   CooksDistance is a measure of scaled change in fitted values. A
%   value larger than three times the mean Cook's distance may be
%   considered influential.
%
%   HatMatrix is an N-by-N matrix H such that Yfit=H*Y where Y is the
%   response vector and Yfit is the vector of fitted response values.
%
%   All of these quantities are computed on the scale of the linear
%   predictor. So, for example, in the equation that defines the hat
%   matrix, 
%
%      Yfit = glm.Fitted.LinearPredictor
%      Y = glm.Fitted.LinearPredictor + glm.Residuals.LinearPredictor
%
%   See also LinearModel, GeneralizedLinearModel, NonLinearModel.
        Diagnostics
    end
    properties(Dependent,GetAccess='protected',SetAccess='protected')
        % "Working" values - created and saved during fit, but subject to
        % being cleared out and recreated when needed
        binomSize_r = []; % reduced, i.e. BinomSize(Subset)
        offset_r = [];    %               Offset(Subset)
    end
    
    methods % get/set methods
        function D = get.Distribution(model)
            name = model.DistributionName;
            [devFun,varFun] = getDistributionFunctions(name);
            D.Name = name;
            D.DevianceFunction = devFun;
            D.VarianceFunction = varFun;
        end
        function link = get.Link(model)
            [linkFun,dlinkFun,ilinkFun] = dfswitchyard('stattestlink',model.Formula.Link,class(model.Coefs));
            if internal.stats.isString(model.Formula.Link)
                link.Name = model.Formula.Link;
            elseif isnumeric(model.Formula.Link);
                link.Name = sprintf('%g',model.Formula.Link);
            else
                link.Name = '';
            end
            link.Link = linkFun;
            link.Derivative = dlinkFun;
            link.Inverse = ilinkFun;
        end
        
        function binomSize_r = get.binomSize_r(model)
            if isempty(model.WorkingValues)
                subset = model.ObservationInfo.Subset;
                model.WorkingValues.binomSize_r = model.ObservationInfo.BinomSize(subset);
            else
                binomSize_r = model.WorkingValues.binomSize_r;
            end
        end
        function offset_r = get.offset_r(model)
            if isempty(model.WorkingValues)
                subset = model.ObservationInfo.Subset;
                offset = model.ObservationInfo.Offset;
                if ~isempty(offset) && ~isscalar(offset)
                    offset = offset(subset);
                end
                model.WorkingValues.offset_r = offset;
            else
                offset_r = model.WorkingValues.offset_r;
            end
        end
        function R = get.Residuals(model)
            R = get_residuals(model);
        end
        
% The following code is removed because it causes a bad interaction with
% the array editor. As a result, the Diagnostics propety does not appear in
% the array editor view of the LinearModel object. Diagnostics property
% access from the command line is provided in the subsref method.

%         function D = get.Diagnostics(model)
%             D = get_diagnostics(model);
%         end
    end % get/set methods
    
    methods(Hidden=true, Access='public')
        function model = GeneralizedLinearModel(varargin)
            if nargin == 0 % special case
                model.Formula = classreg.regr.LinearFormula;
                return
            end
            error(message('stats:GeneralizedLinearModel:NoConstructor'));
        end
        
        % Implementation of VariableEditorPropertyProvider to customize 
        % the display of properties in the Variable Editor
        function isVirtual = isVariableEditorVirtualProp(~,prop)
            % Return true for the Diagnostics property to enable the
            % Variable Editor to derive the Diagnostics property display 
            % without actually accessing the Diagnostics property
            % (which may cause memory overflows).
            isVirtual = strcmp(prop,'Diagnostics');
        end
        function isComplex = isVariableEditorComplexProp(~,~)
            % Diagnostics property should not be complex
            isComplex = false;
        end
        function isSparse = isVariableEditorSparseProp(~,~)
            % Diagnostics property should not be sparse
            isSparse = false;
        end 
        function className = getVariableEditorClassProp(~,~)
            % Diagnostics property in the Variable Editor is table object
            className = 'table';
        end
        function sizeArray = getVariableEditorSize(this,~)
            sizeArray = [size(this.ObservationInfo.Subset,1); 3];
        end      
    end
    
    methods(Access='public')
        function model = step(model,varargin)
%STEP Selectively add or remove terms from a regression model.
%   M2 = STEP(M1) refines the regression model M1 by taking one step of a
%   stepwise regression, and returns the new model as M2. STEP first tries
%   to add a new term that will reduce the AIC value. If none is found, it
%   tries to remove a term if doing so will reduce the AIC value. If none
%   is found, it returns M2 with the same terms as in M1.
%
%   M2 = STEP(M1,'PARAM1',val1,'PARAM2',val2,...) specifies one or more of
%   the following name/value pairs:
%
%      'Lower'     Lower model of terms that must remain in the model,
%                  default='constant'
%      'Upper'     Upper model of terms available to be added to the model,
%                  default='interactions'
%      'Criterion' Criterion to use in evaluating terms to add or remove,
%                  chosen from 'Deviance' (default), 'AIC', 'BIC', 'SSE',
%                  'RSquared', 'AdjRsquared'
%      'PEnter'    For the 'SSE' criterion, a value E such that a term may
%                  be added if its p-value is less than or equal to E. For
%                  the other criteria, a term may be added if the
%                  improvement in the criterion is at least E.
%      'PRemove'   For the 'SSE' criterion, a value R such that a term may
%                  be removed if its p-value is greater or equal to R. For
%                  the other criteria, a term may be added if it reduces
%                  the criterion no more than R.
%      'NSteps'    Maximum number of steps to take, default=1
%      'Verbose'   An integer from 0 to 2 controlling the display of
%                  information. Verbose=1 (the default) displays the action
%                  taken at each step. Verbose=2 also displays the actions
%                  evaluated at each step. Verbose=0 suppresses all
%                  display.
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
%      % Fit a Poisson regression model using random data and a single
%      % predictor, then step other predictors in one at a time.
%      x = randn(100,20);
%      mu = exp(x(:,[5 10 15])*[.4;.2;.3] + 1);
%      y = poissrnd(mu);
%      glm = fitglm(x,y,'y ~ x1','distr','poisson')
%      glm = step(glm)
%      glm = step(glm)
%
%   See also GeneralizedLinearModel, FITGLM, STEPWISEGLM

            % Intercept the criterion before calling parent class method,
            % and translate 'Deviance' to one of two specific types
            [crit,~,args] = internal.stats.parseArgs({'criterion'},{'Deviance'},varargin{:});
            if strncmpi(crit,'deviance',length(crit))
                if model.DispersionEstimated
                    crit = 'deviance_f';
                else
                    crit = 'deviance_chi2';
                end
            end
            model = step@classreg.regr.TermsRegression(model,'Criterion',crit,args{:});
            checkDesignRank(model);
        end
        
        % --------------------------------------------------------------------
        function disp(model)
%DISP Display a GeneralizedLinearModel.
%   DISP(GLM) displays the GeneralizedLinearModel GLM.
%
%   See also GeneralizedLinearModel.
            isLoose = strcmp(get(0,'FormatSpacing'),'loose');
            if (isLoose), fprintf('\n'); end
            fprintf(getString(message('stats:GeneralizedLinearModel:display_GLM')));

            % Fit the formula string to the command window width
            indent = '    ';
            maxWidth = matlab.desktop.commandwindow.size; maxWidth = maxWidth(1) - 1;
            f = model.Formula;
            if any(strcmpi(model.VariableInfo.Class{model.RespLoc},{'nominal','ordinal','categorical'}))
                a = model.VariableInfo.Range(model.RespLoc);
                b = a{:};
                y = char(b(1,2));
                if ~isa(model.Data,'dataset') && ~isa(model.Data,'table')
                    fstr = ['logit(P(Y=''',y,''')) ~ ',model.Formula.LinearPredictor];
                else
                    fstr = ['logit(P(',model.ResponseName,'=''',y,''')) ~ ',model.Formula.LinearPredictor];
                end
            else
                fstr = char(f,maxWidth-length(indent));
            end
            disp([indent fstr]);
            fprintf('%s',getString(message('stats:GeneralizedLinearModel:display_Distribution',indent,model.DistributionName)));
            
            if model.IsFitFromData
                fprintf('%s',getString(message('stats:GeneralizedLinearModel:display_EstimatedCoefficients')));
                disp(model.Coefficients);
                fprintf('%s',getString(message('stats:GeneralizedLinearModel:display_NumObsDFE',model.NumObservations,model.DFE)));
                if model.DispersionEstimated
                    fprintf('%s',getString(message('stats:GeneralizedLinearModel:display_EstimatedDispersion',num2str(model.Dispersion,'%.3g'))));
                else
                    fprintf('%s',getString(message('stats:GeneralizedLinearModel:display_Dispersion',num2str(model.Dispersion,'%.3g'))));
                end
                if hasConstantModelNested(model) && model.NumPredictors > 0
                    d = devianceTest(model);
                    if model.DispersionEstimated
                        fprintf('%s',getString(message('stats:GeneralizedLinearModel:display_FTest',num2str(d.FStat(2),'%.3g'),num2str(d.pValue(2),'%.3g'))));
                    else
                        fprintf('%s',getString(message('stats:GeneralizedLinearModel:display_ChiTest',num2str(d.chi2Stat(2),'%.3g'),num2str(d.pValue(2),'%.3g'))));
                    end
                end
            else
                fprintf('%s',getString(message('stats:GeneralizedLinearModel:display_Coefficients')));
                if any(model.Coefficients.SE > 0)
                    disp(model.Coefficients(:,{'Value' 'SE'}));
                else
                    disp(model.Coefficients(:,{'Value'}));
                end
                if model.Dispersion > 0
                    fprintf('\n%s',getString(message('stats:GeneralizedLinearModel:display_Dispersion',num2str(model.Dispersion,'%.3g'))));
                end
            end
        end
        
        % --------------------------------------------------------------------
        function [ypred,yCI] = predict(model,varargin)
%predict Compute predicted values given predictor values.
%   YPRED = PREDICT(GLM,DS) computes a vector YPRED of predicted values
%   from the GeneralizedLinearModel GLM using predictor variables from the
%   dataset/table DS. DS must contain all of the predictor variables used to
%   create GLM.
%
%   YPRED = PREDICT(GLM,X), where X is a data matrix with the same number
%   of columns as the matrix used to create GLM, computes predictions using
%   the values in X.
%
%   [YPRED,YCI] = PREDICT(...) also returns the two-column matrix YCI
%   containing 95% confidence intervals for the predicted values. These are
%   non-simultaneous intervals for predicting the mean response at the
%   specified predictor values. The lower limits of the bounds are in
%   column 1, and the upper limits are in column 2.
%
%   [...] = PREDICT(GLM,DS,PARAM1,VAL1,PARAM2,VAL2,...) or
%   [...] = PREDICT(GLM,X,PARAM1,VAL1,PARAM2,VAL2,...) specifies one or more
%   of the following name/value pairs: 
%
%      'Alpha'        A value between 0 and 1 to specify the confidence
%                     level as 100(1-ALPHA)%.  Default is 0.05 for 95%
%                     confidence.
%      'Simultaneous' Either true for simultaneous bounds, or false (the
%                     default) for non-simultaneous bounds.
%      'BinomialSize' The value of the binomial N parameter for each
%                     row in DS or X. May be a vector or scalar. The
%                     default value 1 produces YPRED values that are
%                     predicted proportions. This parameter is used only
%                     if GLM is fit to a binomial distribution.
%      'Offset'       Value of the offset for each row in DS or X. May be
%                     a vector or scalar. The offset is used as an
%                     additional predictor with a coefficient value
%                     fixed at 1.0.
%
%   Example:
%      % Fit a logistic model, and compute predicted probabilities and
%      % confidence intervals for the first three observations
%      load hospital
%      formula = 'Smoker ~ Age*Weight*Sex - Age:Weight:Sex';
%      glm = fitglm(hospital,formula,'distr','binomial')
%      [ypred,confint] = predict(glm,hospital(1:3,:))
%
%   See also GeneralizedLinearModel, random.

            if nargin > 1 && ~internal.stats.isString(varargin{1})
                Xpred = varargin{1};
                varargin = varargin(2:end);
                design = designMatrix(model,Xpred);
                offset = 0;
            else
                design = model.Design;
                offset = model.Offset;
            end
            
            paramNames = {'BinomialSize' 'Confidence' 'Simultaneous' 'Alpha' 'Offset'};
            paramDflts = {           []          .95          false  0.05    offset};
            [Npred,conf,simOpt,alpha,offset,supplied] = ...
                internal.stats.parseArgs(paramNames, paramDflts, varargin{:});
            if supplied.Confidence && supplied.Alpha
                error(message('stats:GeneralizedLinearModel:PredictArgCombination', 'Confidence', 'Alpha'))
            end
            if supplied.Alpha
                conf = 1-alpha;
            end
            
            if supplied.BinomialSize && ~isempty(Npred) && ~strcmp(model.DistributionName,'Binomial')
                error(message('stats:GeneralizedLinearModel:NotBinomial'));
            end
            if strcmpi(model.DistributionName,'binomial')
                if supplied.BinomialSize
                    binomSizePV = {'Size' Npred};
                else
                    binomSizePV = {'Size' ones(size(design,1),1)};
                end
            else
                binomSizePV = {};
            end
            
            if nargout < 2
                ypred = glmval(model.Coefs,design,model.Formula.Link,'Constant','off','Offset',offset,binomSizePV{:});
            else
                [R,sigma] = corrcov(model.CoefficientCovariance);
                stats = struct('se',sigma, 'coeffcorr',R, ...
                    'dfe',model.DFE, 's',model.Dispersion, 'estdisp',model.DispersionEstimated);
                [ypred,dylo,dyhi] = ...
                    glmval(model.Coefs,design,model.Formula.Link,stats,'Constant','off', ...
                           'Confidence',conf,binomSizePV{:},'simul',simOpt,'Offset',offset);
                yCI = [ypred-dylo ypred+dyhi];
            end
        end
        
        % --------------------------------------------------------------------
        function ysim = random(model,ds,varargin)
%random Generate random response values given predictor values.
%   YNEW = RANDOM(GLM,DS) generates a vector YNEW of random values from the
%   GeneralizedLinearModel GLM using predictor variables from the
%   dataset/table DS. DS must contain all of the predictor variables used to create GLM.
%   The output YNEW is computed by creating predicted values and generating
%   response values from the fitted distribution.
%
%   YNEW = RANDOM(GLM,X), where X is a data matrix with the same number of
%   columns as the matrix used to create GLM, generates responses using the
%   values in X.
%
%   For binomial and Poisson fits, the YNEW values are generated with the
%   specified distribution. There is no adjustment for any estimated
%   dispersion.
%
%   YNEW = RANDOM(GLM,DS,PARAM1,VAL1,PARAM2,VAL2,...) or
%   YNEW = RANDOM(GLM,X,PARAM1,VAL1,PARAM2,VAL2,...) specifies one or more
%   of the following name/value pairs: 
%
%      'BinomialSize' The value of the binomial N parameter for each
%                     row in DS or X. May be a vector or scalar. The
%                     default value is 1. This parameter is used only
%                     if GLM is fit to a binomial distribution.
%      'Offset'       Value of the offset for each row in DS or X. May be
%                     a vector or scalar. The offset is used as an
%                     additional predictor with a coefficient value
%                     fixed at 1.0.
%
%   Example:
%      % Generate some Poisson data
%      x = 2 + randn(100,1);
%      mu = exp(1 + x/2);
%      y = poissrnd(mu);
%      subplot(1,2,1)
%      scatter(x,y)
%
%      % Fit a regression model and generate new random response values
%      % from that model
%      glm = fitglm(x,y,'y ~ x1','distr','poisson')
%      y2 = random(glm,x);
%      subplot(1,2,2)
%      scatter(x,y2)
%
%   See also GeneralizedLinearModel, predict, Dispersion.

            if strcmpi(model.DistributionName,'binomial')
                % Compute predicted probabilities, use N parameter later
                paramNames = {'BinomialSize'};
                paramDflts = {           1  };
                [N,supplied,otherargs] = ...
                internal.stats.parseArgs(paramNames, paramDflts, varargin{:});
                varargin = otherargs;
                if supplied.BinomialSize && ~isequal(N,1)
                    varargin(end+1:end+2) = {'BinomialSize' 1};
                end
                if isempty(N)
                    N = 1;
                elseif ~isscalar(N) && numel(N)~=size(ds,1)
                    error(message('stats:GeneralizedLinearModel:RandomNSize', size(ds,1)))
                end
            end
            ypred = predict(model,ds,varargin{:});
            switch lower(model.DistributionName)
            case 'binomial'
                if any(strcmpi(model.VariableInfo.Class{model.RespLoc},{'nominal','ordinal','categorical'}))
                    % Categorical responses
                    if N~=1
                        error(message('stats:GeneralizedLinearModel:BadBinomialSize'))
                    else                        
                        a = model.VariableInfo.Range(model.RespLoc);
                        b = a{:};
                        y = binornd(N,ypred);
                        y1(y==0) = b(1,1);
                        y1(y==1) = b(1,2);
                        y1 = y1';
                        ysim=y1;
                    end
                else
                    ysim = binornd(N,ypred);
                end
                
            case 'gamma'
                ysim = gamrnd(1./model.Dispersion,ypred.*model.Dispersion);
            case 'inverse gaussian'
                ysim = random('inversegaussian',ypred,1./model.Dispersion);
            case 'normal'
                ysim = normrnd(ypred,sqrt(model.Dispersion));
            case 'poisson'
                ysim = poissrnd(ypred);
            end
        end

        % --------------------------------------------------------------------
        function tbl = devianceTest(model)
%devianceTest Analysis of deviance
%   TBL = devianceTest(GLM) displays an analysis of deviance table for the
%   GeneralizedLinearModel GLM. The table displays a test of whether the
%   fitted model fits significantly better than a constant model.
%
%   Example:
%      % Generate some Poisson data
%      x = 2 + randn(100,1);
%      mu = exp(1 + x/2);
%      y = poissrnd(mu);
%
%      % Fit a regression model and examine the significance of the fit
%      glm = fitglm(x,y,'y ~ x1','distr','poisson')
%      devianceTest(glm)
%
%   See also GeneralizedLinearModel.
            dev = [model.DevianceNull; model.Deviance];
            dfe = [model.NumObservations-1; model.DFE];
            dfr = model.NumEstimatedCoefficients - 1;
            
            if model.DispersionEstimated
                statName = 'FStat';
                stat = -diff(dev) ./ (dfr * model.Dispersion); % F statistic
                p = fcdf(1./stat,dfe(2),dfr); % fpval(fstat,dfr,dfe);
            else
                statName = 'chi2Stat';
                stat = -diff(dev); % chi-sqared statistic
                p = gammainc(stat/2,dfr/2,'upper'); % chi2pval(chi2stat,dfr);
            end
            
            if ~hasConstantModelNested(model)
                warning(message('stats:GeneralizedLinearModel:NoIntercept'));
            end
            
            f0 = model.Formula; f0.Terms = zeros(1,model.NumVariables);
            tbl = table(dev, ...
                          dfe, ...
                          internal.stats.DoubleTableColumn([NaN; stat],[true; false]), ...
                          internal.stats.DoubleTableColumn([NaN; p],[true; false]), ...
                          'VariableNames',{'Deviance' 'DFE' statName 'pValue'}, ...
                          'RowNames',{char(f0,40) char(model.Formula,40)});
            tbl.Properties.DimensionNames = {'Fits' 'Variables'};
        end
        
    function hout = plotDiagnostics(model,plottype,varargin)
%plotDiagnostics Plot diagnostics of fitted model
%    plotDiagnostics(GLM,PLOTTYPE) plots diagnostics from
%    GeneralizedLinearModel GLM in a plot of type PLOTTYPE. The default
%    value of PLOTTYPE is 'leverage'. Valid values for PLOTTYPE are:
%
%       'contour'      residual vs. leverage with overlaid Cook's contours
%       'cookd'        Cook's distance
%       'leverage'     leverage (diagonal of Hat matrix)
%
%    H = plotDiagnostics(...) returns handles to the lines in the plot.
%
%    The PLOTTYPE argument can be followed by parameter/value pairs to
%    specify additional properties of the primary line in the plot. For
%    example, plotDiagnostics(GLM,'cookd','Marker','s') uses a square
%    marker.
%
%    The data cursor tool in the figure window will display the X and Y
%    values for any data point, along with the observation name or number.
%    It also displays the coefficient name for 'dfbetas'.
%
%    See also NonLinearModel, plotResiduals.
            if nargin<2
                plottype = 'leverage';
            else
                alltypes = {'contour' 'cookd' 'leverage'};
                plottype = internal.stats.getParamVal(plottype,alltypes,'PLOTTYPE');
            end
            f = classreg.regr.modelutils.plotDiagnostics(model,plottype,varargin{:});
            if nargout>0
                hout = f;
            end
        end
        
    function hout = plotResiduals(model,plottype,varargin)
%plotResiduals Plot residuals of fitted model
%    plotResiduals(MODEL,PLOTTYPE) plots the residuals from model MODEL in
%    a plot of type PLOTTYPE. Valid values for PLOTTYPE are:
%
%       'caseorder'     residuals vs. case (row) order
%       'fitted'        residuals vs. fitted values
%       'histogram'     histogram (default)
%       'lagged'        residual vs. lagged residual (r(t) vs. r(t-1))
%       'probability'   normal probability plot
%       'symmetry'      symmetry plot
%
%    plotResiduals(MODEL,PLOTTYPE,'ResidualType',RTYPE) plots the residuals
%    of type RTYPE, which can be any of the following:
%
%       'Raw'             Observed minus fitted values
%       'LinearPredictor' Residuals on the linear predictor scale, equal to
%                         the adjusted response value minus the fitted
%                         linear combination of the predictors.
%       'Pearson'         Residuals divided by the estimated standard
%                         deviation of the response.
%       'Anscombe'        Residuals defined on transformed data with the
%                         transformation chosen to remove skewness.
%       'Deviance'        Residuals based on the contribution of each
%                         observation to the deviance.
%
%    H = plotResiduals(...) returns a handle to the lines or patches in the
%    plot.
%
%    The PLOTTYPE argument can be followed by parameter/value pairs to
%    specify additional properties of the primary line in the plot. For
%    example, plotResiduals(M,'fitted','Marker','s') uses a square marker.
%
%    For many of these plots, the data cursor tool in the figure window
%    will display the X and Y values for any data point, along with the
%    observation name or number.
%
%    See also GeneralizedLinearModel.

            if nargin<2
                plottype = 'histogram';
            end
            [residtype,~,args] = internal.stats.parseArgs({'residualtype'},{'Raw'},varargin{:});
            varargin = args;
            residtype = internal.stats.getParamVal(residtype,...
                         GeneralizedLinearModel.SupportedResidualTypes,'''ResidualType''');

            f = classreg.regr.modelutils.plotResiduals(model,plottype,'ResidualType',residtype,varargin{:});
            if nargout>0
                hout = f;
            end
        end
    
    function fout = plotSlice(model)
%plotSlice Plot slices through the fitted regression surface.
%   plotSlice(GLM) creates a new figure containing a series of plots, each
%   representing a slice through the regression surface predicted by GLM.
%   For each plot, the surface slice is shown as a function of a single
%   predictor variable, with the other predictor variables held constant.
%
%   If the model GLM includes an offset term, that term is set to 0 in
%   computing the regression surface.
%
%   See also GeneralizedLinearModel, predict.

            f = classreg.regr.modelutils.plotSlice(model);
            
            % This model doesn't support bounds for a new observation, so
            % disable the menu items that select that
            h1 = findobj(f, 'Tag','boundsCurve');
            h2 = findobj(f, 'Tag','boundsObservation');
            set([h1 h2],'Enable','off');
            
            if nargout>0
                fout = f;
            end
        end
    end % public
    
    methods(Hidden=true)
        % --------------------------------------------------------------------
        function t = title(model)
            strLHS = linkstr(model.Formula.Link,model.ResponseName);
            strFunArgs = internal.stats.strCollapse(model.Formula.PredictorNames,',');
            t = sprintf('%s = glm(%s)',strLHS,strFunArgs);
        end
        % --------------------------------------------------------------------
        function v = varianceParam(model)
            v = model.Dispersion;
        end
    end        

    methods(Access='protected')
        function tbl = tstats(model)
            if model.DispersionEstimated
                nobs = model.NumObservations;
            else
                nobs = Inf;
            end
            tbl = classreg.regr.modelutils.tstats(model.Coefs,sqrt(diag(model.CoefficientCovariance)), ...
                                      nobs,model.CoefficientNames);
        end

        function model = assignData(model,X,y,w,offset,binomN,asCat,dummyCoding,varNames,excl)
            model = assignData@classreg.regr.TermsRegression(model,X,y,w,asCat,dummyCoding,varNames,excl);
            y = getResponse(model);
            n = length(y);
            
            % Response is ordinarily a vector
            if strcmpi(model.DistributionName,'binomial')
                % Binomial response may be a vector y with n=1 assumed,
                % or a vector y with args   'BinomialSize',n,
                % or a matrix [y,n]
                if isvector(y)
                    if isempty(binomN)
                        binomN = ones(n,1);
                    elseif ~isnumeric(binomN) || ~all(binomN==round(binomN)) ...
                                          || any(binomN<0) ...
                                          || ~(isscalar(binomN) ...
                                              || (isvector(binomN) && numel(binomN)==n))
                        error(message('stats:GeneralizedLinearModel:BadNValues'));
                    elseif isscalar(binomN)
                        binomN = repmat(binomN,n,1);
                    end
                    model.ObservationInfo.BinomSize = binomN(:);
                elseif size(y,2)==2 && isnumeric(y) ...
                                    && all(all(binomN==round(binomN))) ...
                                    && all(all(binomN>=0))
                    model.ObservationInfo.BinomSize = y(:,2);
                    model.ResponseType = 'BinomialCounts';
                else
                    error(message('stats:GeneralizedLinearModel:BadBinomialResponse'));
                end
            elseif ~isvector(y)
                error(message('stats:GeneralizedLinearModel:BadResponse'));
            end
            if isempty(offset)
                offset = zeros(n,1);
            elseif ~isnumeric(offset) || ~(isscalar(offset) || ...
                                           (isvector(offset) && numel(offset)==n))
                error(message('stats:GeneralizedLinearModel:BadOffset', n))
            elseif isscalar(offset)
                offset = repmat(offset,n,1);
            end
            model.Offset = offset(:);
        end
        
        % --------------------------------------------------------------------
        function model = selectObservations(model,exclude,missing)
            if nargin < 3, missing = []; end
            model = selectObservations@classreg.regr.TermsRegression(model,exclude,missing);
            subset = model.ObservationInfo.Subset;
            if strcmpi(model.DistributionName,'binomial')
                % Populate fields in the WorkingValues structure
                model.WorkingValues.binomSize_r = model.ObservationInfo.BinomSize(subset);
            end
            offset = model.Offset;
            if ~isempty(offset) && ~isscalar(offset)
                offset = offset(subset);
            end
            model.WorkingValues.offset_r = offset;
        end
        
        % --------------------------------------------------------------------
        function model = fitter(model)
            X = getData(model);
            [model.Design,model.CoefTerm,model.CoefficientNames] = designMatrix(model,X);
            
            % Populate the design_r field in the WorkingValues structure
            model.WorkingValues.design_r = create_design_r(model);
            
            if strcmpi(model.DistributionName,'binomial')
                response_r = [model.y_r model.binomSize_r];
            else
                response_r = model.y_r;
            end
            estDisp = model.DispersionEstimated || ...
                      ~any(strcmpi(model.DistributionName,{'binomial' 'poisson'}));
            [model.Coefs,model.Deviance,stats] = ...
                glmfit(model.design_r,response_r,model.DistributionName,'EstDisp',estDisp, ...
                       'Constant','off','Link',model.Formula.Link,'Weights',model.w_r,...
                       'Offset',model.offset_r,'RankWarn',false);
            model.CoefficientCovariance = stats.covb;
            model.DFE = stats.dfe;
            model.Dispersion = stats.s^2;

            subset = model.ObservationInfo.Subset;
            wts = NaN(size(subset));
            wts(subset) = stats.wts;
            model.IRLSWeights = wts;
        end
        
        % --------------------------------------------------------------------
        function model = postFit(model)
            % Do housework after fitting
            model = postFit@classreg.regr.TermsRegression(model);
            if strcmpi(model.DistributionName,'binomial')
                y = [model.y_r model.binomSize_r];
            else
                y = model.y_r;
            end
            estDisp = model.DispersionEstimated && ...
                      ~any(strcmpi(model.DistributionName,{'binomial' 'poisson'}));
            [~,model.DevianceNull] = ...
                glmfit(ones(size(y,1),1),y,model.DistributionName,'EstDisp',estDisp, ...
                       'Constant','off','Link',model.Formula.Link,'Weights',model.w_r,...
                       'Offset',model.offset_r);

        end
        
%         % --------------------------------------------------------------------
%         function binomSize_r = create_binomSize_r(model)
%             binomSize_r = model.ObservationInfo.BinomSize(model.ObservationInfo.Subset);
%         end
%         
        % --------------------------------------------------------------------
        function ypred = predictPredictorMatrix(model,Xpred)
            % Assume the columns correspond in order to the predictor variables
            design = designMatrix(model,Xpred,true);
            ypred = glmval(model.Coefs,design,model.Formula.Link,'Constant','off');
        end
        
        % --------------------------------------------------------------------
        function L = getlogLikelihood(model)
            subset = model.ObservationInfo.Subset;
            switch lower(model.DistributionName)
            case 'binomial'
                pHat_r = predict(model); pHat_r = pHat_r(subset);
                L = sum(model.w_r .* binologpdf(model.y_r,model.binomSize_r,pHat_r));
            case 'gamma'
                aHat = 1 ./ model.Dispersion;
                bHat_r = predict(model) .* model.Dispersion; bHat_r = bHat_r(subset);
                L = sum(model.w_r .* gamlogpdf(model.y_r,aHat,bHat_r));
            case 'inverse gaussian'
                muHat_r = predict(model); muHat_r = muHat_r(subset);
                lambdaHat = 1 ./ model.Dispersion;
                L = sum(model.w_r .* invglogpdf(model.y_r,muHat_r,lambdaHat));
            case 'normal'
                muHat_r = predict(model); muHat_r = muHat_r(subset);
                sigmaHat = sqrt(model.DFE/model.NumObservations*model.Dispersion);
                L = sum(model.w_r .* normlogpdf(model.y_r,muHat_r,sigmaHat));
            case 'poisson'
                lambdaHat_r = predict(model); lambdaHat_r = lambdaHat_r(subset);
                L = sum(model.w_r .* poisslogpdf(model.y_r,lambdaHat_r));
            end
        end
        
        % --------------------------------------------------------------------
        function L0 = logLikelihoodNull(model)
            % We want normalized weights for calculating weighted means
            % (The var() function normalizes internally.)
            w_r_normalized = model.w_r / sum(model.w_r);
            switch lower(model.DistributionName)
            case 'binomial'
                p0 = sum(model.w_r .* model.y_r) ./ sum(model.w_r .* model.binomSize_r);
                L0 = sum(model.w_r .* binologpdf(model.y_r,model.binomSize_r,p0));
            case 'gamma'
                mu0 = sum(w_r_normalized .* model.y_r);
                b0 = var(model.y_r,model.w_r) ./ mu0;
                a0 = mu0 ./ b0;
                L0 = sum(model.w_r .* gamlogpdf(model.y_r,a0,b0));
            case 'inverse gaussian'
                mu0 = sum(w_r_normalized .* model.y_r);
                lambda0 = (mu0^3) ./ var(model.y_r,model.w_r); 
                L0 = sum(model.w_r .* invglogpdf(model.y_r,mu0,lambda0));
            case 'normal'
                mu0 = sum(w_r_normalized .* model.y_r);
                sigma0 = std(model.y_r,model.w_r);
                L0 = sum(model.w_r .* normlogpdf(model.y_r,mu0,sigma0));
            case 'poisson'
                lambda0 = sum(w_r_normalized .* model.y_r);
                L0 = sum(model.w_r .* poisslogpdf(model.y_r,lambda0));
            end
        end
        
        % --------------------------------------------------------------------
        function yfit = get_fitted(model,type)
            if nargin < 2
                Response = get_fitted(model,'response');
                LinearPredictor = get_fitted(model,'linearpredictor');
                yfit = table(Response,LinearPredictor,'RowNames',model.ObservationNames);
                if strcmpi(model.DistributionName,'binomial')
                    yfit.Probability = get_fitted(model,'probability');
                end
            else
                yfit = predict(model);
                switch lower(type)
                case 'response'
                    if strcmpi(model.DistributionName,'binomial')
                        yfit = yfit .* model.ObservationInfo.BinomSize;
                    end
                case 'probability'
                    if ~strcmpi(model.DistributionName,'binomial')
                        error(message('stats:GeneralizedLinearModel:NoProbabilityFit'));
                    end
                case 'linearpredictor'
                    link = dfswitchyard('stattestlink',model.Formula.Link,class(yfit));
                    yfit = link(yfit);
                otherwise
                    error(message('stats:GeneralizedLinearModel:UnrecognizedFit', type));
                end
            end
        end
        
        % --------------------------------------------------------------------
        function d = get_CooksDistance(model)
            if hasData(model)
                w = get_CombinedWeights_r(model,false);
                r = get_residuals(model,'linearpredictor');
                h = model.Leverage;
                d = w .* abs(r).^2 .* (h./(1-h).^2)./(model.NumEstimatedCoefficients*varianceParam(model));
            else
                d = [];
            end
        end
        function r = get_residuals(model,type)
            if nargin < 2 % return all residual types in a table array
                Raw = get_residuals(model,'raw');
                LinearPredictor = get_residuals(model,'linearpredictor');
                Pearson = get_residuals(model,'pearson');
                Anscombe = get_residuals(model,'anscombe');
                Deviance = get_residuals(model,'deviance'); %#ok<PROP>
                r = table(Raw,LinearPredictor,Pearson,Anscombe,Deviance, ...
                            'RowNames',model.ObservationNames); %#ok<PROP>
            else % return a single type of residual
                y = getResponse(model);
                if strcmpi(model.DistributionName,'binomial')
                    N = model.ObservationInfo.BinomSize;
                    y = y(:,1) ./ N;
                else
                    N = 1;
                end
                mu = predict(model);
                
                switch lower(type)
                case 'raw'
                    if strcmpi(model.DistributionName,'binomial')
                        r = (y - mu) .* N;
                    else
                        r = y - mu;
                    end
                case 'linearpredictor'
                    dlink = model.Link.Derivative;
                    deta = dlink(mu);
                    % From the glm fitting iterations:
                    %    eta = offset + x * b;
                    %    mu = link(eta);
                    %    z = eta + (y - mu) .* deta;
                    % On the z scale, then, the residual is:
                    r = (y - mu) .* deta;
                case 'pearson'
                    if strcmpi(model.DistributionName,'binomial')
                        r = (y - mu) ./ (sqrt(model.Distribution.VarianceFunction(mu,N)) + (y==mu));
                    else
                        r = (y - mu) ./ (sqrt(model.Distribution.VarianceFunction(mu)) + (y==mu));
                    end
                case 'anscombe'
                    switch(lower(model.DistributionName))
                    case 'normal'
                        r = y - mu;
                    case 'binomial'
                        t = 2/3;
                        r = beta(t,t) * (betainc(y,t,t)-betainc(mu,t,t)) ...
                            ./ ((mu.*(1-mu)).^(1/6) ./ sqrt(N));
                    case 'poisson'
                        r = 1.5 * ((y.^(2/3) - mu.^(2/3)) ./ mu.^(1/6));
                    case 'gamma'
                        r = 3 * (y.^(1/3) - mu.^(1/3)) ./ mu.^(1/3);
                    case 'inverse gaussian'
                        r = (log(y) - log(mu)) ./ mu;
                    end
                case 'deviance'
                    devFun = model.Distribution.DevianceFunction;
                    if strcmpi(model.DistributionName,'binomial')
                        r = sign(y - mu) .* sqrt(max(0,devFun(mu,y,N)));
                    else
                        r = sign(y - mu) .* sqrt(max(0,devFun(mu,y)));
                    end
                otherwise
                    error(message('stats:GeneralizedLinearModel:UnrecognizedResidual', type));
                end
            end
        end
        function w = get_CombinedWeights_r(model,reduce)
            w = model.ObservationInfo.Weights .* model.IRLSWeights;
            if nargin<2 || reduce
                subset = model.ObservationInfo.Subset;
                w = w(subset);
            end
        end
        
        % -------------------- pass throughs to modelutils -------------------
        function crit = get_rsquared(model,type)
            stats = struct('SSE',model.SSE, ...
                           'SST',model.SST, ...
                           'DFE',model.DFE, ...
                           'NumObservations',model.NumObservations, ...
                           'LogLikelihood',model.LogLikelihood, ...
                           'LogLikelihoodNull',model.LogLikelihoodNull, ...
                           'Deviance',model.Deviance, ...
                           'DevianceNull',model.DevianceNull);
            if nargin < 2
                crit = classreg.regr.modelutils.rsquared(stats,'all',true);
            else
                crit = classreg.regr.modelutils.rsquared(stats,type);
            end
        end
    end % protected
        
    methods(Static, Access='public', Hidden)
        function model = fit(X,varargin) % [X, y | DS], modelDef, ...
% Not intended to be called directly. Use FITGLM to fit a GeneralizedLinearModel.
%
%   See also FITGLM.

            [X,y,haveDataset,otherArgs] = GeneralizedLinearModel.handleDataArgs(X,varargin{:});
            
            % VarNames are optional names for the X matrix and y vector.  A
            % dataset/table defines its own list of names, so this is not accepted
            % with a dataset/table.
            
            % PredictorVars is an optional list of the subset of variables to
            % actually use as predictors in the model, and is only needed with
            % an alias.  A terms matrix or a formula string already defines
            % which variables to use without that.  ResponseVar is an optional
            % name that is not needed with a formula string.
            
            paramNames = {'Distribution' 'Link' 'Intercept' 'PredictorVars' 'ResponseVar' ...
                          'Weights' 'Exclude' 'CategoricalVars' 'VarNames' 'DummyVarCoding' ...
                          'BinomialSize' 'Offset' 'DispersionFlag' 'rankwarn'};
            paramDflts = {'normal' '' true [] [] [] [] [] [] 'reference' [] 0 [] true};

            % Default model is main effects only.
            if isempty(otherArgs)
                modelDef = 'linear';
            else
                arg1 = otherArgs{1};
                if mod(length(otherArgs),2)==1 % odd, model followed by pairs
                    modelDef = arg1;
                    otherArgs(1) = [];
                elseif internal.stats.isString(arg1) && ...
                       any(strncmpi(arg1,paramNames,length(arg1)))
                    % omitted model but included name/value pairs
                    modelDef = 'linear';
                end
            end

            [distr,link,intercept,predictorVars,responseVar, ...
             weights,exclude,asCatVar,varNames,dummyCoding,binomN,offset,...
             disperFlag,rankwarn,supplied] = ...
                internal.stats.parseArgs(paramNames, paramDflts, otherArgs{:});
            
            alldist = {'Normal' 'Poisson' 'Binomial', 'Gamma', 'Inverse Gaussian'};
            distr = internal.stats.getParamVal(distr,alldist,'''Distribution''');
            if isempty(disperFlag)
                disperFlag = ~any(strcmpi(distr,{'Binomial', 'Poisson'}));
            end
            
            if ~isempty(binomN) && ~strcmp(distr,'Binomial')
                error(message('stats:GeneralizedLinearModel:NotBinomial'));
            end
            
            model = GeneralizedLinearModel();
            model.DistributionName = distr;
            clink = canonicalLink(distr);
            
            if isscalar(disperFlag) && (islogical(disperFlag) || ...
                    (isnumeric(disperFlag) && (disperFlag==0 || disperFlag==1)))
                model.DispersionEstimated = disperFlag==1;
            else
                error(message('stats:GeneralizedLinearModel:BadDispersion'));
            end
            
   % Categorical responses
            modelType = {'linear';'interactions';'purequadratic';'quadratic';'polyIJK...'};            
            if haveDataset == 1 && isa(X.(X.Properties.VariableNames{end}),'categorical') && any(strcmpi(modelDef,modelType))
                % dataset/table with last column as response
                dsVarNames = X.Properties.VariableNames{end};
                y0 = X.(dsVarNames);
                [X.(X.Properties.VariableNames{end}), classname] = categ2num(y0);
            else 
                y0=[];
            end
            
            model.Formula = GeneralizedLinearModel.createFormula(supplied,modelDef,X, ...
                predictorVars,responseVar,intercept,link,varNames,haveDataset,clink);
            
            if isempty(y0)
                if haveDataset == 0 && isa(y,'categorical')
                    % data in separate X and y
                    y0 = y;
                    [y, classname] = categ2num(y);
                elseif haveDataset == 1 && isa(X.(model.Formula.ResponseName),'categorical')
                    % dataset/table with formula/ResponseName
                    y0 = X.(model.Formula.ResponseName);
                    [X.(model.Formula.ResponseName), classname]= categ2num(y0);
                else
                    y0 = 0;
                end
            end
            
            model = assignData(model,X,y,weights,offset,binomN,asCatVar,dummyCoding,model.Formula.VariableNames,exclude);
            silent = classreg.regr.LinearFormula.isModelAlias(modelDef);
            model = removeCategoricalPowers(model,silent);
            model = doFit(model);
            
            if isa(y0,'categorical')
                a = class(y0);
                model.VariableInfo.Class{model.RespLoc} = a;
                if strcmpi(model.VariableInfo.Class{model.RespLoc},{'nominal'})
                    model.VariableInfo.Range{model.RespLoc} = nominal(classname');
                elseif strcmpi(model.VariableInfo.Class{model.RespLoc},{'ordinal'})
                    model.VariableInfo.Range{model.RespLoc} = ordinal(classname');
                elseif strcmpi(model.VariableInfo.Class{model.RespLoc},{'categorical'})
                    model.VariableInfo.Range{model.RespLoc} = categorical(classname');
                end
            end
            
            model = updateVarRange(model); % omit excluded points from range
           
            if rankwarn
                checkDesignRank(model);
            end
        end
        
        function model = stepwise(X,varargin) % [X, y | DS], start, ...
% Not intended to be called directly. Use STEPWISEGLM to fit a
% GenerarlizedLinearModel by stepwise regression.
%
%   See also STEPWISEGLM.

            [X,y,haveDataset,otherArgs] = GeneralizedLinearModel.handleDataArgs(X,varargin{:});
            
            if isempty(otherArgs)
                start = 'constant';
            else
                start = otherArgs{1};
                otherArgs(1) = [];
            end
            
            paramNames = {'Distribution' 'Link' 'Intercept' 'PredictorVars' 'ResponseVar' 'Weights' ...
                'Exclude' 'CategoricalVars' 'VarNames' 'Lower' 'Upper' 'Criterion' ...
                'PEnter' 'PRemove' 'NSteps' 'Verbose' 'BinomialSize' 'Offset' 'DispersionFlag'};
            paramDflts = {'normal' [] true [] [] [] [] [] [] 'constant' 'interactions' 'Deviance' [] [] Inf 1};
            [distr,link,intercept,predictorVars,responseVar,weights,exclude, ...
                asCatVar,varNames,lower,upper,crit,penter,premove,nsteps,verbose,binomN,offset,...
                disperFlag,supplied] = ...
                internal.stats.parseArgs(paramNames, paramDflts, otherArgs{:});
            
            if ~isscalar(verbose) || ~ismember(verbose,0:2)
                error(message('stats:LinearModel:BadVerbose'));
            end
            
            alldist = {'Normal' 'Poisson' 'Binomial', 'Gamma', 'Inverse Gaussian'};
            distr = internal.stats.getParamVal(distr,alldist,'''Distribution''');
            if isempty(disperFlag)
                disperFlag = ~any(strcmpi(distr,{'Binomial', 'Poisson'}));
            end
            if ~(isscalar(disperFlag) && (islogical(disperFlag) || ...
                    (isnumeric(disperFlag) && (disperFlag==0 || disperFlag==1))))
                error(message('stats:GeneralizedLinearModel:BadDispersion'));
            end
            
            clink = canonicalLink(distr);
            
            % *** need more reconciliation among start, lower, upper, and between those
            % and link, intercept, predictorVars, varNames
            if ~supplied.ResponseVar && (classreg.regr.LinearFormula.isTermsMatrix(start) || classreg.regr.LinearFormula.isModelAlias(start))
                if isa(lower,'classreg.regr.LinearFormula')
                    responseVar = lower.ResponseName;
                    supplied.ResponseVar = true;
                else
                    if internal.stats.isString(lower) && ~classreg.regr.LinearFormula.isModelAlias(lower)
                        lower = GeneralizedLinearModel.createFormula(supplied,lower,X, ...
                            predictorVars,responseVar,intercept,link,varNames,haveDataset,clink);
                        responseVar = lower.ResponseName;
                        supplied.ResponseVar = true;
                    elseif isa(upper,'classreg.regr.LinearFormula')
                        responseVar = upper.ResponseName;
                        supplied.ResponseVar = true;
                    else
                        if internal.stats.isString(upper) && ~classreg.regr.LinearFormula.isModelAlias(upper)
                            upper = GeneralizedLinearModel.createFormula(supplied,upper,X, ...
                                predictorVars,responseVar,intercept,link,varNames,haveDataset,clink);
                            responseVar = upper.ResponseName;
                            supplied.ResponseVar = true;
                        end
                    end
                end
            end
            
            if ~isa(start,'classreg.regr.LinearFormula')
                start = GeneralizedLinearModel.createFormula(supplied,start,X, ...
                    predictorVars,responseVar,intercept,link,varNames,haveDataset,clink);
            end
            if ~isa(lower,'classreg.regr.LinearFormula')
                if classreg.regr.LinearFormula.isModelAlias(lower)
                    if supplied.PredictorVars
                        lower = {lower,predictorVars};
                    end
                end
                lower = classreg.regr.LinearFormula(lower,start.VariableNames,start.ResponseName,start.HasIntercept,start.Link);
            end
            if ~isa(upper,'classreg.regr.LinearFormula')
                if classreg.regr.LinearFormula.isModelAlias(upper)
                    if supplied.PredictorVars
                        upper = {upper,predictorVars};
                    end
                end
                upper = classreg.regr.LinearFormula(upper,start.VariableNames,start.ResponseName,start.HasIntercept,start.Link);
            end
            
            if haveDataset
                model = GeneralizedLinearModel.fit(X,start.Terms,'Distribution',distr, ...
                    'Link',start.Link,'ResponseVar',start.ResponseName,'Weights',weights, ...
                    'Exclude',exclude,'CategoricalVars',asCatVar,...
                    'Disper',disperFlag,'BinomialSize',binomN,'Offset',offset,'RankWarn',false);
            else
                model = GeneralizedLinearModel.fit(X,y,start.Terms,'Distribution',distr, ...
                    'Link',start.Link,'ResponseVar',start.ResponseName,'Weights',weights, ...
                    'Exclude',exclude,'CategoricalVars',asCatVar,'VarNames',start.VariableNames,...
                    'Disper',disperFlag,'BinomialSize',binomN,'Offset',offset,'RankWarn',false);
            end
            
            model.Steps.Start = start;
            model.Steps.Lower = lower;
            model.Steps.Upper = upper;
            if strcmpi(crit,'deviance')
                if model.DispersionEstimated
                    model.Steps.Criterion = 'deviance_f';
                else
                    model.Steps.Criterion = 'deviance_chi2';
                end
            else
                model.Steps.Criterion = crit;
            end
            model.Steps.PEnter = penter;
            model.Steps.PRemove = premove;
            model.Steps.History = [];
            
            model = stepwiseFitter(model,nsteps,verbose);
            checkDesignRank(model);
        end
    end % static hidden
    
    methods(Static, Access='protected')
        function [addTest,addThreshold,removeTest,removeThreshold,reportedNames,testName] = getStepwiseTests(crit,Steps)
            localCrits = {'deviance_f' 'deviance_chi2'};
            if internal.stats.isString(crit) && ismember(crit,localCrits)
                addThreshold = Steps.PEnter;
                removeThreshold = Steps.PRemove;
                if isempty(addThreshold), addThreshold = 0.05; end
                if isempty(removeThreshold), removeThreshold = 0.10; end
                
                if addThreshold>=removeThreshold
                    error(message('stats:LinearModel:BadSmallerThreshold', sprintf( '%g', addThreshold ), sprintf( '%g', removeThreshold ), crit));
                end

                switch lower(crit)
                    case 'deviance_f'
                        addTest = @(proposed,current) f_test(proposed,current,'up');
                        removeTest = @(proposed,current) f_test(current,proposed,'down');
                        reportedNames = {'Deviance' 'FStat' 'PValue'};
                        testName = 'pValue';
                    case 'deviance_chi2'
                        addTest = @(proposed,current) chi2_test(proposed,current,'up');
                        removeTest = @(proposed,current) chi2_test(current,proposed,'down');
                        reportedNames = {'Deviance' 'Chi2Stat' 'PValue'};
                        testName = 'pValue';
                end
            else
                [addTest,addThreshold,removeTest,removeThreshold,reportedNames,testName] ...
                    = classreg.regr.TermsRegression.getStepwiseTests(crit,Steps);
            end
        end % static protected
    end
end

function [p,pReported,reportedVals] = f_test(fit1,fit0,direction)
dev1 = fit1.Deviance; % the larger model
dfDenom = fit1.DFE;
if isempty(fit0)
    dev0 = dev1;
    dfNumer = 0;
else
    dev0 = fit0.Deviance; % the smaller model
    dfNumer = fit0.DFE - fit1.DFE;
end
F = ((dev0-dev1)/dfNumer) / fit1.Dispersion;
p = fcdf(1./F,dfDenom,dfNumer); % fpval(F,dfNumer,dfDenom)
pReported = p;
if strcmp(direction,'up') % testing significance of larger model
    reportedVals = {dev1 F p};
else % 'down', testing non-significance of larger model
    reportedVals = {dev0 F p};
end
end

function [p,pReported,reportedVals] = chi2_test(fit1,fit0,direction)
dev1 = fit1.Deviance; % the larger model
if isempty(fit0)
    dev0 = dev1;
    df = 0;
else
    dev0 = fit0.Deviance; % the smaller model
    df = fit0.DFE - fit1.DFE;
end
x2 = (dev0-dev1) / df;
p = gammainc(x2/2,df/2,'upper'); % chi2pval(x2,df)
pReported = p;
if strcmp(direction,'up') % testing significance of larger model
    reportedVals = {dev1 x2 p};
else % 'down', testing non-significance of larger model
    reportedVals = {dev0 x2 p};
end
end

%-------------------------------------------------------------------------
function str = linkstr(link,responseName)
switch lower(link)
    case 'identity',   str = sprintf('%s',responseName);
    case {'inverse' 'reciprocal'}, str = sprintf('1/%s',responseName);
    case 'log',        str = sprintf('log(%s)',responseName);
    case 'logit',      str = sprintf('log(%s/(1-%s))',responseName,responseName);
    case 'probit',     str = sprintf('probit(%s)',responseName);
    case 'loglog',     str = sprintf('log(-log(%s))',responseName);
    case 'comploglog', str = sprintf('log(-log(1-%s))',responseName);
    otherwise
        if isnumeric(link)
            str = sprintf('(%s)^%d',responseName,link);
        else
            str = sprintf('CustomLink(%s)',responseName);
        end
end
end


% ----------------------------------------------------------------------------
function logy = binologpdf(x,n,p)
logy = gammaln(n+1) - gammaln(x+1) - gammaln(n-x+1) + x.*log(p) + (n-x).*log1p(-p);
end
function logy = gamlogpdf(x,a,b)
z = x./b;
logy = (a-1).*log(z) - z - gammaln(a) - log(b);
end
function logy = invglogpdf(x,mu,lambda)
logy = .5*log(lambda./(2.*pi.*x.^3)) + (-0.5.*lambda.*(x./mu - 2 + mu./x)./mu);
end
function logy = normlogpdf(x,mu,sigma)
logy = (-0.5 * ((x - mu)./sigma).^2) - log(sqrt(2*pi) .* sigma);
end
function logy = poisslogpdf(x,lambda)
logy = (-lambda + x.*log(lambda) - gammaln(x+1));
end

function [devFun,varFun] = getDistributionFunctions(name)
switch lower(name)
    case 'normal'
        devFun = @(mu,y) (y - mu).^2;
        varFun = @(mu) ones(size(mu),class(mu));
    case 'binomial'
        devFun = @(mu,y,N) 2*N.*(y.*log((y+(y==0))./mu) + (1-y).*log((1-y+(y==1))./(1-mu)));
        varFun = @(mu,N) mu.*(1-mu) ./ N;
    case 'poisson'
        devFun = @(mu,y) 2*(y .* (log((y+(y==0)) ./ mu)) - (y - mu));
        varFun = @(mu) mu;
    case 'gamma'
        devFun = @(mu,y) 2*(-log(y ./ mu) + (y - mu) ./ mu);
        varFun = @(mu) mu.^2;
    case 'inverse gaussian'
        devFun = @(mu,y) ((y - mu)./mu).^2 ./  y;
        varFun = @(mu) mu.^3;
end
end

% --------------------------------------------------------------------
function link = canonicalLink(distr)
switch lower(distr)
    case 'normal'
        link = 'identity';
    case 'binomial'
        link = 'logit';
    case 'poisson'
        link = 'log';
    case 'gamma'
        link = 'reciprocal';
    case 'inverse gaussian'
        link = -2;
end
end

% --------------------------------------------------------------------
function [y, classname] = categ2num(x)  % change two-level categorical response to 0 and 1
[y, classname] = grp2idx(x);
nc = length(classname);
if nc > 2
    error(message('stats:glmfit:TwoLevelCategory'));
end
y(y==1) = 0;
y(y==2) = 1;
end
