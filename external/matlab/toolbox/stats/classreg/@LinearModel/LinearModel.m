classdef (Sealed = true) LinearModel < classreg.regr.TermsRegression
%LinearModel Fitted multiple linear regression model.
%   LM = FITLM(...) fits a linear model to data. The fitted model LM is a
%   LinearModel that can predict a response as a linear function of
%   predictor variables and terms created from predictor variables.
%
%   LinearModel methods:
%       addTerms - Add terms to linear model
%       removeTerms - Remove terms from linear model
%       step - Selectively add or remove terms from linear model
%       anova - Analysis of variance
%       coefCI - Coefficient confidence intervals
%       coefTest - Linear hypothesis test on coefficients
%       predict - Compute predicted values given predictor values
%       feval - Evaluate model as a function
%       random - Generate random response values given predictor values
%       dwtest - Durbin-Watson test for autocorrelation in residuals
%       plot - Summary plot of regression model
%       plotAdded - Plot of marginal effect of a single term
%       plotAdjustedResponse - Plot of response and one predictor
%       plotDiagnostics - Plot of regression diagnostics
%       plotEffects - Plot of main effects of predictors
%       plotInteraction - Plot of interaction effects of two predictors
%       plotResiduals - Plot of residuals
%       plotSlice - Plot of slices through fitted regression surface
%
%   LinearModel properties:
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
%       Diagnostics - Regression diagnostics
%       MSE - Mean squared error (estimate of residual variance)
%       RMSE - Root mean squared error (estimate of residual standard deviation)
%       DFE - Degrees of freedom for residuals
%       SSE - Error sum of squares
%       SST - Total sum of squares
%       SSR - Regression sum of squares
%       Steps - Stepwise fit results
%       Robust - Robust fit results
%       Formula - Representation of the model used in this fit
%       LogLikelihood - Log of likelihood function at coefficient estimates
%       CoefficientCovariance - Covariance matrix for coefficient estimates
%       CoefficientNames - Coefficient names
%       NumCoefficients - Number of coefficients
%       NumEstimatedCoefficients - Number of estimated coefficients
%
%   See also FITLM, GeneralizedLinearModel, NonLinearModel, STEPWISELM.

%   Copyright 2011-2013 The MathWorks, Inc.

    properties(Constant,Hidden)
        SupportedResidualTypes = {'Raw' 'Pearson' 'Standardized' 'Studentized'};
    end
    properties(GetAccess='public',SetAccess='protected')
        
%MSE - Mean squared error.
%   The MSE property is the mean squared error. It is an estimate of the
%   variance of the error term.
%
%   See also LinearModel, anova.
        MSE
        
%Robust - Robust regression results.
%   If the model was fit using robust regression, the Robust property
%   is a structure containing information about that fit. If the model was
%   not fit using robust regression, this property is empty.
%
%   The Robust structure contains the following fields:
%      RobustWgtFun  Weight function used for the fit, default 'bisquare'
%      Tune          Tuning constant used
%      Weights       Robust weights used at the final iteration
%
%   See also LinearModel.
        Robust = [];
    end
    properties(GetAccess='protected',SetAccess='protected')
        Qy
        R
        Rtol
        Q
    end
    properties(Dependent=true,GetAccess='public',SetAccess='protected')
%Residuals - Residual values.
%   The Residuals property is a table of residuals. It is a table array that
%   has one row for each observation and the following variables:
%
%        'Raw'          Observed minus fitted values
%        'Pearson'      Raw residuals divided by RMSE
%        'Standardized' Raw residuals divided by their estimated standard
%                       deviation
%        'Studentized'  Raw residuals divided by an independent (delete-1)
%                       estimate of their standard deviation
%
%   To obtain any of these columns as a vector, index into the property
%   using dot notation. For example, in the model M, the ordinary or
%   raw residual vector is
%
%      r = M.Residuals.Raw
%
%   See also LinearModel, plotResiduals, Fitted, predict, random.
        Residuals

%Fitted - Fitted (predicted) values.
%   The Fitted property is a vector of fitted values.
%
%   The fitted values are computed using the predictor values used to fit
%   the model. Use the PREDICT method to compute predictions for other
%   predictor values and to compute confidence bounds for the predicted
%   values.
%
%   See also LinearModel, Residuals, predict, random.
        Fitted

%Diagnostics - Regression diagnostics.
%   The Diagnostics property is a structure containing a set of diagnostics
%   helpful in finding outliers and influential observations. Many describe
%   the effect on the fit of deleting single observations.  The structure
%   contains the following fields:
%      Leverage  Diagonal elements of the Hat matrix
%      Dffits    Scaled change in fitted values with row deletion
%      CooksDistance Cook's measure of scaled change in fitted values
%      S2_i      Residual variance estimate with row deletion
%      Dfbetas   Scaled change in coefficient estimates with row deletion
%      CovRatio  Covariance determinant ratio with row deletion
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
%   Dffits is the scaled change in the fitted values for each observation
%   that would result from excluding that observation from the fit. Values
%   with an absolute value larger than 2*sqrt(P/N) may be considered
%   influential.
%
%   CooksDistance is another measure of scaled change in fitted values. A
%   value larger than three times the mean Cook's distance may be
%   considered influential.
%
%   S2_i is a set of residual variance estimates obtained by deleting each
%   observation in turn. These can be compared with the value of the MSE
%   property.
%
%   Dfbetas is an N-by-P matrix of the scaled change in the coefficient
%   estimates that would result from excluding each observation in turn.
%   Values larger than 3/sqrt(N) in absolute value indicate that the
%   observation has a large influence on the corresponding coefficient.
%
%   CovRatio is the ratio of the determinant of the coefficient covariance
%   matrix with each observation deleted in turn to the determinant of the
%   covariance matrix for the full model. Values larger than 1+3*P/N or
%   smaller than 1-3*P/N indicate influential points.
%
%   HatMatrix is an N-by-N matrix H such that Yfit=H*Y where Y is the
%   response vector and Yfit is the vector of fitted response values.
%
%   See also LinearModel, GeneralizedLinearModel, NonLinearModel.
        Diagnostics

%RMSE - Root mean squared error.
%   The RMSE property is the root mean squared error. It is an estimate of
%   the standard deviation of the error term.
%
%   See also anova.
        RMSE
    end
    
    methods % get/set methods
        function yfit = get.Fitted(model)
            yfit = get_fitted(model);
        end
        function r = get.Residuals(model)
            r = get_residuals(model);
        end
        function s = get.RMSE(model)
            s = sqrt(model.MSE);
        end
        
% The following code is removed because it causes a bad interaction with
% the array editor. As a result, the Diagnostics propety does not appear in
% the array editor view of the LinearModel object. Diagnostics property
% access from the command line is provided in the subsref method.

        function D = get.Diagnostics(model)
            D = get_diagnostics(model);
        end
    end
    methods(Access='protected')
        function s2_i = get_S2_i(model)
            r = getResponse(model) - predict(model);
            h = model.Leverage; % from parent class
            wt = get_CombinedWeights_r(model,false);
            delta_i = wt .* abs(r).^2 ./ (1-h);
            if any(h==1)
                % If any points are completely determined by their own
                % observation, then removing those points doesn't decrease
                % the SSE and also doesn't decrease the DFE
                newdf = repmat(model.DFE-1,length(h),1);
                delta_i(h==1) = 0;
                newdf(h==1) = newdf(h==1) + 1;
            else
                newdf = model.DFE-1;
            end
            s2_i = max(0,model.SSE - delta_i) ./ newdf;
            subset = model.ObservationInfo.Subset;
            s2_i(~subset & ~isnan(s2_i)) = 0;
        end
        function dfbetas = get_Dfbetas(model)
            rows = model.ObservationInfo.Subset;
            w_r = get_CombinedWeights_r(model);
            [~,~,~,R1,~,~,~,Q1] = lsfit(model.design_r,model.y_r,w_r);
            C = Q1/R1';
            e_i = model.Residuals.Studentized(rows,:);
            h = model.Leverage(rows,:); % from parent class
            dfbetas = zeros(length(e_i),size(C,2));
            dfb = bsxfun(@rdivide,C,sqrt(sum(C.^2)));
            dfb = bsxfun(@times,dfb, sqrt(w_r).*e_i./sqrt(1-h));
            dfbetas(rows,:) = dfb;
        end
        function dffits = get_Dffits(model)
            e_i = model.Residuals.Studentized;
            wt = get_CombinedWeights_r(model,false);
            h = model.Leverage; % from parent class
            dffits = sqrt(h./(1-h)).*sqrt(wt).*e_i;
        end
        function covr = get_CovRatio(model)
            n = model.NumObservations;
            p = model.NumEstimatedCoefficients;
            wt = get_CombinedWeights_r(model,false);
            e_i = model.Residuals.Studentized;
            h = model.Leverage; % from parent class
            covr = 1 ./ ((((n-p-1+wt.*abs(e_i).^2)./(n-p)).^p).*(1-h));
        end
        function w = get_CombinedWeights_r(model,reduce)
            w = model.ObservationInfo.Weights;
            if ~isempty(model.Robust)
                w = w .* model.Robust.Weights;
            end
            if nargin<2 || reduce
                subset = model.ObservationInfo.Subset;
                w = w(subset);
            end
        end
    end % get/set methods
    
    methods(Hidden=true, Access='public')
        function model = LinearModel(varargin) % modelDef, coefs, ...
            if nargin == 0 % special case
                model.Formula = classreg.regr.LinearFormula;
                return
            end
            error(message('stats:LinearModel:NoConstructor'));
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
            sizeArray = [size(this.ObservationInfo.Subset,1); 7];
        end      
    end
    
    methods(Access='public')
        function disp(model)
%DISP Display a LinearModel.
%   DISP(LM) displays the LinearModel LM.
%
%   See also LinearModel.
            isLoose = strcmp(get(0,'FormatSpacing'),'loose');
            if (isLoose), fprintf('\n'); end
            if isempty(model.Robust)
                fprintf('%s',getString(message('stats:LinearModel:display_LinearRegressionModel')));
            else
                fprintf('%s',getString(message('stats:LinearModel:display_LinearRegressionModelrobustFit')));
            end

            % Fit the formula string to the command window width
            indent = '    ';
            maxWidth = matlab.desktop.commandwindow.size; maxWidth = maxWidth(1) - 1;
            f = model.Formula;
            fstr = char(f,maxWidth-length(indent));
            disp([indent fstr]);
            
            if model.IsFitFromData
                fprintf('%s',getString(message('stats:LinearModel:display_EstimatedCoefficients')));
                disp(model.Coefficients);
                fprintf('%s',getString(message('stats:LinearModel:display_NumObservationsDFE',model.NumObservations,model.DFE)));
                fprintf('%s',getString(message('stats:LinearModel:display_RMSE',num2str(model.RMSE,'%.3g'))));
                if hasConstantModelNested(model) && model.NumPredictors > 0
                    rsq = get_rsquared(model,{'ordinary','adjusted'});
                    fprintf('%s',getString(message('stats:LinearModel:display_RsquaredAdj',num2str(rsq(1),'%.3g'),num2str(rsq(2),'%.3g'))));
                    [f,p] = fTest(model);
                    fprintf('%s',getString(message('stats:LinearModel:display_Ftest',num2str(f,'%.3g'),num2str(p,'%.3g'))));
                end
            else
                fprintf(getString(message('stats:LinearModel:display_Coefficients')));
                if any(model.CoefSE > 0)
                    disp(model.Coefficients(:,{'Value' 'SE'}));
                else
                    disp(model.Coefficients(:,{'Value'}));
                end
                if model.MSE > 0
                    fprintf('\n%s',getString(message('stats:LinearModel:display_RMSE',num2str(model.RMSE,'%.3g'))));
                end
            end
        end
        
        % --------------------------------------------------------------------
        function model = step(model,varargin)
%STEP Selectively add or remove terms from a regression model.
%   M2 = STEP(M1) refines the regression model M1 by taking one step of a
%   stepwise regression, and returns the new model as M2. STEP first tries
%   to add a new term that will reduce the AIC value. If none is found, it
%   tries to remove a term if doing so will reduce the AIC value. If none
%   is found, it returns M2 with the same terms as in M1.
%
%   The STEP method is not available with robust fits.
%
%   M2 = STEP(M1,'PARAM1',val1,'PARAM2',val2,...) specifies one or more of
%   the following name/value pairs:
%
%      'Lower'     Lower model of terms that must remain in the model,
%                  default='constant'
%      'Upper'     Upper model of terms available to be added to the model,
%                  default='interactions'
%      'Criterion' Criterion to use in evaluating terms to add or remove,
%                  chosen from 'SSE' (default) 'AIC', 'BIC', 'RSquared',
%                  'AdjRsquared'
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
%      'SSE'         0.05   < 0.10       p-value for F test
%      'AIC'         0      < 0.01       change in AIC
%      'BIC'         0      < 0.01       change in BIC
%      'Rsquared'    0.1    > 0.05       increase in R-squared
%      'AdjRsquared' 0      > -0.05      increase in adjusted R-squared
%
%    Example:
%       % Fit model to car data; check for any term to add or remove from a
%       % quadratic model
%       load carsmall
%       d = dataset(MPG,Weight);
%       d.Year = ordinal(Model_Year);
%       lm1 = fitlm(d,'MPG ~ Year + Weight + Weight^2')
%       lm2 = step(lm1,'upper','quadratic')
%
%   See also LinearModel, FITLM, STEPWISELM.
            if ~isempty(model.Robust)
                error(message('stats:LinearModel:NoRobustStepwise'));
            end
            model = step@classreg.regr.TermsRegression(model,varargin{:});
            checkDesignRank(model);
        end
        
        % --------------------------------------------------------------------
        function [ypred,yCI] = predict(model,varargin)
%predict Compute predicted values given predictor values.
%   YPRED = PREDICT(LM,DS) computes a vector YPRED of predicted values from
%   the LinearModel LM using predictor variables from the dataset/table DS. DS
%   must contain all of the predictor variables used to create LM.
%
%   YPRED = PREDICT(LM,X), where X is a data matrix with the same number of
%   columns as the matrix used to create LM, computes predictions using the
%   values in X.
%
%   [YPRED,YCI] = PREDICT(...) also returns the two-column matrix YCI
%   containing 95% confidence intervals for the predicted values. These are
%   non-simultaneous intervals for predicting the mean response at the
%   specified predictor values. The lower limits of the bounds are in
%   column 1, and the upper limits are in column 2.
%
%   [...] = PREDICT(LM,DS,PARAM1,VAL1,PARAM2,VAL2,...) or
%   [...] = PREDICT(LM,X,PARAM1,VAL1,PARAM2,VAL2,...) specifies one or more
%   of the following name/value pairs: 
%
%      'Alpha'        A value between 0 and 1 to specify the confidence
%                     level as 100(1-ALPHA)%.  Default is 0.05 for 95%
%                     confidence.
%      'Simultaneous' Either true for simultaneous bounds, or false (the
%                     default) for non-simultaneous bounds.
%      'Prediction'   Either 'curve' (the default) to compute confidence
%                     intervals for the curve (function value) at X, or
%                     'observation' for prediction intervals for a new
%                     observation at X.
%
%   Example:
%      % Create a regression model and use it to compute predictions
%      % and confidence intervals for the value of the function for
%      % the first three observations
%      load carsmall
%      d = dataset(MPG,Weight);
%      d.Year = ordinal(Model_Year);
%      lm = fitlm(d,'MPG ~ Year + Weight + Weight^2')
%      [fitted,confint] = predict(lm,d(1:3,:))
%
%   See also LinearModel, random.

            if nargin > 1 && ~internal.stats.isString(varargin{1})
                Xpred = varargin{1};
                varargin = varargin(2:end);
                design = designMatrix(model,Xpred);
            else
                design = model.Design;
            end
            
            paramNames = {'Alpha' 'Simultaneous' 'Prediction' 'Confidence'};
            paramDflts = {    .05          false      'curve' .95 };
            [alpha,simOpt,predOpt,conf,supplied] = ...
                internal.stats.parseArgs(paramNames, paramDflts, varargin{:});
            if supplied.Confidence && supplied.Alpha
                error(message('stats:LinearModel:ArgCombination'))
            end
            if supplied.Confidence
                alpha = 1-conf;
            end

            predOpt = internal.stats.getParamVal(predOpt,...
                           {'curve' 'observation'},'''Prediction''');
            predOpt = strcmpi(predOpt,'observation');
            simOpt = internal.stats.parseOnOff(simOpt,'''Simultaneous''');
            if nargout < 2
                ypred = predci(design,model.Coefs);
            else
                [ypred, yCI] = predci(design,model.Coefs,...
                                  model.CoefficientCovariance,model.MSE,...
                                  model.DFE,alpha,simOpt,predOpt,model.Formula.HasIntercept);
            end
        end
        
        % --------------------------------------------------------------------
        function ysim = random(model,varargin) % model, Xsim, Nsim
%random Generate random response values given predictor values.
%   YNEW = RANDOM(LM,DS) generates a vector YNEW of random values from
%   the LinearModel LM using predictor variables from the dataset/table DS. 
%   DS must contain all of the predictor variables used to create LM. The
%   output YNEW is computed by creating predicted values and adding new
%   random noise values with standard deviation LM.RMSE.
%
%   YNEW = RANDOM(LM,X), where X is a data matrix with the same number of
%   columns as the matrix used to create LM, generates responses using the
%   values in X.
%
%   Example:
%      % Plot car MPG as a function of Weight
%      load carsmall
%      subplot(1,2,1);
%      scatter(Weight,MPG)
% 
%      % Same plot using new random responses from a fitted model
%      subplot(1,2,2);
%      lm = fitlm(Weight,MPG,'quadratic');
%      mrand = random(lm,Weight)
%      scatter(Weight,mrand)
%
%   See also LinearModel, predict, RMSE.
            ypred = predict(model,varargin{:});
            ysim = normrnd(ypred,model.RMSE);
        end
        
        % -------------------- pass throughs to modelutils -------------------
function hout = plotDiagnostics(model,varargin)
%plotDiagnostics Plot diagnostics of fitted model
%    plotDiagnostics(LM,PLOTTYPE) plots diagnostics from LinearModel LM in
%    a plot of type PLOTTYPE. The default value of PLOTTYPE is 'dffits'.
%    Valid values for PLOTTYPE are:
%
%       'contour'      residual vs. leverage with overlayed Cook's contours
%       'cookd'        Cook's distance
%       'covratio'     delete-1 ratio of determinant of covariance
%       'dfbetas'      scaled delete-1 coefficient estimates
%       'dffits'       scaled delete-1 fitted values
%       'leverage'     leverage (diagonal of Hat matrix)
%       's2_i'         delete-1 variance estimate
%
%    H = plotDiagnostics(...) returns handles to the lines in the plot.
%
%    The PLOTTYPE argument can be followed by parameter/value pairs to
%    specify additional properties of the primary line in the plot. For
%    example, plotDiagnostics(LM,'cookd','Marker','s') uses a square
%    marker.
%
%    The data cursor tool in the figure window will display the X and Y
%    values for any data point, along with the observation name or number.
%    It also displays the coefficient name for 'dfbetas'.
%
%    Example:
%      % Plot the leverage in a fitted regression model
%      load carsmall
%      d = dataset(MPG,Weight);
%      d.Year = ordinal(Model_Year);
%      lm = fitlm(d,'MPG ~ Year + Weight + Weight^2')
%      plotDiagnostics(lm,'leverage')
%
%      % Look at the data for the high-leverage points, and note that
%      % their Weight values are near the extremes
%      high = find(lm.Diagnostics.Leverage>0.11)
%      d(high,:)
%
%    See also LinearModel, plotResiduals.
            f = classreg.regr.modelutils.plotDiagnostics(model,varargin{:});
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
%        'Raw'          Observed minus fitted values
%        'Pearson'      Raw residuals divided by RMSE
%        'Standardized' Raw residuals divided by their estimated standard
%                       deviation
%        'Studentized'  Raw residuals divided by an independent (delete-1)
%                       estimate of their standard deviation
%
%    H = plotResiduals(...) returns a handle to the lines or patches in the
%    plot.
%
%    The PLOTTYPE or RTYPE arguments can be followed by parameter/value
%    pairs to specify additional properties of the primary line in the
%    plot. For example, plotResiduals(M,'fitted','Marker','s') uses a
%    square marker.
%
%    For many of these plots, the data cursor tool in the figure window
%    will display the X and Y values for any data point, along with the
%    observation name or number.
%
%    Example:
%      % Make a normal probability plot of the raw residuals in a fitted
%      % regression model
%      load carsmall
%      d = dataset(MPG,Weight);
%      d.Year = ordinal(Model_Year);
%      lm = fitlm(d,'MPG ~ Year + Weight + Weight^2')
%      plotResiduals(lm,'probability')
%
%      % Examine the points with the largest residuals
%      high = find(lm.Residuals.Raw > 8)
%      d(high,:)
%
%    See also LinearModel, plotDiagnostics.
            if nargin<2
                plottype = 'histogram';
            end
            [residtype,~,args] = internal.stats.parseArgs({'residualtype'},{'Raw'},varargin{:});
            varargin = args;
            residtype = internal.stats.getParamVal(residtype,...
                         LinearModel.SupportedResidualTypes,'''ResidualType''');
            internal.stats.plotargchk(varargin{:});

            f = classreg.regr.modelutils.plotResiduals(model,plottype,'ResidualType',residtype,varargin{:});
            if nargout>0
                hout = f;
            end
        end
    
    function fout = plotSlice(model)
%plotSlice Plot slices through the fitted regression surface.
%   plotSlice(LM) creates a new figure containing a series of plots, each
%   representing a slice through the regression surface predicted by LM.
%   For each plot, the surface slice is shown as a function of a single
%   predictor variable, with the other predictor variables held constant.
%
%   If there are more than eight predictors, plotSlice selects the first
%   five of them for plotting. You can use the Predictors menu to control
%   which predictors are plotted.
%
%    Example:
%      % Make a slice plot for a fitted regression model
%      load carsmall
%      d = dataset(MPG,Weight);
%      d.Year = ordinal(Model_Year);
%      lm = fitlm(d,'MPG ~ Year + Weight + Weight^2')
%      plotSlice(lm)
%
%   See also LinearModel, predict, RMSE.
            f = classreg.regr.modelutils.plotSlice(model);
            if nargout>0
                fout = f;
            end
        end
        
        % ------------------- Externally-defined methods ---------------------
        [p,dw] = dwtest(model,option,tail)
    end % public
    
    methods(Access='protected')
        function model = fitter(model)
            X = getData(model);
            [model.Design,model.CoefTerm,model.CoefficientNames] = designMatrix(model,X);
            
            % Populate the design_r field in the WorkingValues structure
            model.WorkingValues.design_r = create_design_r(model);
            
            if isempty(model.Robust)
                [model.Coefs,model.MSE,model.CoefficientCovariance,model.R,model.Qy,model.DFE,model.Rtol] ...
                    = lsfit(model.design_r,model.y_r,model.w_r);
            else
                [model.Coefs,stats] ...
                    = robustfit(model.design_r,model.y_r,model.Robust.RobustWgtFun,model.Robust.Tune,'off',model.w_r,false);
                model.CoefficientCovariance = stats.covb;
                model.DFE = stats.dfe;
                model.MSE = stats.s^2;
                model.Rtol = stats.Rtol;

                w = NaN(size(model.ObservationInfo,1),1);
                w(model.ObservationInfo.Subset) = stats.w;
                model.Robust.Weights = w;

                model.R = stats.R;
                model.Qy = stats.Qy;
            end
        end
        function model = postFit(model)
            % Do housework after fitting
            model = postFit@classreg.regr.TermsRegression(model);
            
            % Override SSE and SST to take any robust fitting into account
            model.SSE = model.DFE * model.MSE;
            model.SST = model.SSR + model.SSE;
        end
        
        % --------------------------------------------------------------------
        function ypred = predictPredictorMatrix(model,Xpred)
            % Assume the columns correspond in order to the predictor variables
            ypred = designMatrix(model,Xpred,true) * model.Coefs;
        end
        
        % --------------------------------------------------------------------
        function L = getlogLikelihood(model)
            subset = model.ObservationInfo.Subset;
            muHat_r = predict(model);
            muHat_r = muHat_r(subset);
            sigmaHat = sqrt(model.DFE/model.NumObservations*model.MSE);
            if sigmaHat==0
                % We have a perfect fit, so an infinite log likelihood
                L = Inf;
            else
                L = sum(model.w_r .* normlogpdf(model.y_r,muHat_r,sigmaHat));
            end
        end
        
        % --------------------------------------------------------------------
        function L0 = logLikelihoodNull(model)
            mu0 = sum(model.w_r .* model.y_r) / sum(model.w_r);
            sigma0 = std(model.y_r,model.w_r);
            L0 = sum(model.w_r .* normlogpdf(model.y_r,mu0,sigma0));
        end
        
        % --------------------------------------------------------------------
        function D = get_diagnostics(model,type)
            if nargin<2 % return all diagnostics in a table
                HatMatrix = get_diagnostics(model,'hatmatrix');
                CooksDistance = get_diagnostics(model,'cooksdistance');
                Dffits = get_diagnostics(model,'dffits');
                S2_i = get_diagnostics(model,'s2_i');
                Dfbetas = get_diagnostics(model,'dfbetas');
                CovRatio = get_diagnostics(model,'covratio');
                Leverage = model.Leverage;
                D = table(Leverage,CooksDistance,...
                            Dffits,S2_i,CovRatio,Dfbetas,HatMatrix,...
                            'RowNames',model.ObservationNames);
            else        % return a single diagnostic
                subset = model.ObservationInfo.Subset;
                switch(lower(type))
                    case 'leverage'
                        D = model.Leverage;
                        D(~subset,:) = 0;                  
                    case 'hatmatrix'
                        try
                            D = get_HatMatrix(model);
                        catch ME
                            warning(message('stats:LinearModel:HatMatrixError', ...
                                ME.message));
                            D = zeros(length(subset),0);
                        end
                        D(~subset,:) = 0;
                    case 'cooksdistance'
                        D = get_CooksDistance(model);
                        D(~subset,:) = NaN;
                    case 'dffits'
                        D = get_Dffits(model);
                        D(~subset,:) = NaN;
                    case 's2_i'
                        D = get_S2_i(model);
                        D(~subset,:) = NaN;
                    case 'dfbetas'
                        D = get_Dfbetas(model);
                        D(~subset,:) = 0;                  
                    case 'covratio'
                        D = get_CovRatio(model);
                        D(~subset,:) = NaN;
                    otherwise
                        error(message('stats:LinearModel:UnrecognizedDiagnostic', type));
                end
            end
        end
        function r = get_residuals(model,type)
            if nargin < 2 % return all residual types in a table array
                Raw = get_residuals(model,'raw');
                Pearson = get_residuals(model,'pearson');
                Studentized = get_residuals(model,'studentized');
                Standardized = get_residuals(model,'standardized');
                r = table(Raw,Pearson,Studentized,Standardized, ...
                          'RowNames',model.ObservationNames);
            else % return a single type of residual
                subset = model.ObservationInfo.Subset;
                raw = getResponse(model) - predict(model);
                switch lower(type)
                case 'raw'
                    r = raw;
                case 'pearson'
                    r = raw ./ model.RMSE;
                case 'studentized' % "externally studentized", using Delete 1 Variances
                    h = model.Leverage; % from parent class
                    s2_i = get_S2_i(model);
                    r = raw ./ sqrt(s2_i .* (1-h));
                case 'standardized' % "internally studentized", using MSE
                    h = model.Leverage; % from parent class
                    r = raw ./ (sqrt(model.MSE * (1-h)));
                otherwise
                    error(message('stats:LinearModel:UnrecognizedResidual', type));
                end
                r(~subset) = NaN;
            end
        end
    end % protected
    methods(Hidden=true, Access='public') % public to allow testing
        [fxi,fxiVar] = getAdjustedResponse(model,var,xi,terminfo)
        [effects,effectSEs,effectnames,effectXs] = getEffects(model,vars,terminfo)
        [effect,effectSE,effectName] = getConditionalEffect(model,var1,var2,xi1,terminfo)

        % --------------------------------------------------------------------
        function t = title(model)
            strLHS = model.ResponseName;
            strFunArgs = internal.stats.strCollapse(model.Formula.PredictorNames,',');
            t = sprintf( '%s = lm(%s)',strLHS,strFunArgs);
        end
        % --------------------------------------------------------------------
        function v = varianceParam(model)
            v = model.MSE;
        end
    end % hidden public

    methods(Static, Access='public', Hidden)
        function model = fit(X,varargin)
% Not intended to be called directly. Use FITLM to fit a LinearModel.
%
%   See also FITLM.

            [X,y,haveDataset,otherArgs] = LinearModel.handleDataArgs(X,varargin{:});
            
            % VarNames are optional names for the X matrix and y vector.  A
            % dataset/table defines its own list of names, so this is not accepted
            % with a dataset/table.
            
            % PredictorVars is an optional list of the subset of variables to
            % actually use as predictors in the model, and is only needed with
            % an alias.  A terms matrix or a formula string already defines
            % which variables to use without that.  ResponseVar is an optional
            % name that is not needed with a formula string.
            
            % rankwarn is undocumented and is used during stepwise fitting
            
            paramNames = {'Intercept' 'PredictorVars' 'ResponseVar' ...
                          'Weights' 'Exclude' 'CategoricalVars' 'VarNames'...
                          'RobustOpts' 'DummyVarCoding' 'rankwarn'};
            paramDflts = {[] [] [] [] [] [] [] [] 'reference' true};

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
            
            [intercept,predictorVars,responseVar,weights,exclude, ...
             asCatVar,varNames,robustOpts,dummyCoding,rankwarn,supplied] = ...
                internal.stats.parseArgs(paramNames, paramDflts, otherArgs{:});
            
            model = LinearModel();
            
            model.Robust = classreg.regr.FitObject.checkRobust(robustOpts);
            model.Formula = LinearModel.createFormula(supplied,modelDef,X, ...
                predictorVars,responseVar,intercept,varNames,haveDataset);
            model = assignData(model,X,y,weights,asCatVar,dummyCoding,model.Formula.VariableNames,exclude);
            
            silent = classreg.regr.LinearFormula.isModelAlias(modelDef);
            model = removeCategoricalPowers(model,silent);
            
            model = doFit(model);
            
            model = updateVarRange(model); % omit excluded points from range

            if rankwarn
                checkDesignRank(model);
            end
        end
        
       
        % ------------------- Externally-defined methods ---------------------
        model = stepwise(X,varargin)
    end % static public hidden
    
    methods(Static, Access='protected')
        function formula = createFormula(supplied,modelDef,X,predictorVars,responseVar,intercept,varNames,haveDataset)
            supplied.Link = false;
            formula = classreg.regr.TermsRegression.createFormula(supplied,modelDef, ...
                X,predictorVars,responseVar,intercept,'identity',varNames,haveDataset);
        end
    end % static protected
end

% ----------------------------------------------------------------------------
function logy = normlogpdf(x,mu,sigma)
logy = (-0.5 * ((x - mu)./sigma).^2) - log(sqrt(2*pi) .* sigma);
end


% ----------------------------------------------------------------------------
function [ypred, yci] = predci(X,beta,Sigma,mse,dfe,alpha,sim,pred,hasintercept)

% Compute the predicted values at the new X.
ypred = X * beta;

if nargout > 1 % Calculate confidence interval

    if (pred) % prediction interval for new observations
        varpred = sum((X*Sigma) .* X,2) + mse;
    else % confi interval for fitted curve
        varpred = sum((X*Sigma) .* X,2);
    end

    if (sim) % simultaneous
        if (pred)
            % For new observations.
            if (hasintercept)
                % Jacobian has constant column.
                sch = length(beta);
            else
                % Need to use a conservative setting.
                sch = length(beta) + 1;
            end
        else
            % For fitted curve.
            sch = length(beta);
        end
       crit = sqrt(sch * finv(1-alpha, sch, dfe));
    else % pointwise
       crit = tinv(1-alpha/2,dfe);
    end
    delta = sqrt(varpred) * crit;
    yci = [ypred-delta ypred+delta];
end
end

% ----------------------------------------------------------------------------
function [b,mse,S,R1,Qy1,dfe,Rtol,Q1] = lsfit(X,y,w)
% LSFIT Weighted least squares fit

%   Copyright 2011 The MathWorks, Inc.

[nobs,nvar] = size(X); % num observations, num predictor variables

% Weights not given, assume equal.
if nargin < 3 || isempty(w)
    w = [];

% Weights given.
elseif isvector(w) && numel(w)==nobs && all(w>=0)
    D = sqrt(w(:));
    X = bsxfun(@times,D,X);
    y = bsxfun(@times,D,y);
    % w is OK

else
    error(message('stats:LinearModel:InvalidWeights', nobs));
end

outClass = superiorfloat(X,y,w);

% Factor the design matrix and transform the response vector.
[Q,R,perm] = qr(X,0);
Qy = Q'*y;

% Use the rank-revealing QR to remove dependent columns of X.
if isempty(R)
    Rtol = 1;
    keepCols = zeros(1,0);
else
    Rtol = abs(R(1)).*max(nobs,nvar).*eps(class(R));
    if isrow(R)
        keepCols = 1;
    else
        keepCols = find(abs(diag(R)) > Rtol);
    end
end

rankX = length(keepCols);
R0 = R;
perm0 = perm;
if rankX < nvar
    R = R(keepCols,keepCols);
    Qy = Qy(keepCols,:);
    perm = perm(keepCols);
end

% Compute the LS coefficients, filling in zeros in elements corresponding
% to rows of R that were thrown out.
b = zeros(nvar,1,outClass);
b(perm,1) = R \ Qy;

if nargout > 1
    % Compute the MSE.
    dfe = nobs - rankX;
    if dfe > 0
        sst = sum(y.*conj(y),1);
        ssx = sum(Qy.*conj(Qy),1);
        mse = max(0, sst-ssx) ./ dfe;
    else % rankX == nobs, and so Xb == y exactly
        mse = zeros(1,1,outClass);
    end

    % Compute the covariance matrix of the LS estimates.  Fill in zeros
    % corresponding to exact zero coefficients.
    Rinv = R \ eye(rankX,outClass);
    if nargout > 2
        S = zeros(nvar,nvar,outClass);
        S(perm,perm) = Rinv*Rinv' .* mse; % guaranteed to be hermitian
    end
    
    % Return unpermuted, unreduced versions of Q*y and R
    if nargout > 3
        Qy1 = zeros(nvar,1);
        Qy1(perm,1) = Qy;
        R1 = zeros(nvar,nvar,outClass);
        R1(perm,perm0) = R0(keepCols,:);
        Q1 = zeros(size(X),outClass);
        Q1(:,perm) = Q(:,keepCols);
    end
end
end
