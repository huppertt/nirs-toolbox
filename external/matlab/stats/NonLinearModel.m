classdef (Sealed = true) NonLinearModel < classreg.regr.ParametricRegression
%NonLinearModel Fitted nonlinear regression model.
%   NLM = FITNLM(...) fits a nonlinear function to data. The fitted model
%   NLM is a NonLinearModel that can predict a response as a nonlinear
%   function of predictor variables.
%
%   NonLinearModel methods:
%       coefCI - Coefficient confidence intervals
%       coefTest - Linear hypothesis test on coefficients
%       predict - Compute predicted values given predictor values
%       feval - Evaluate model as a function
%       random - Generate random response values given predictor values
%       plotDiagnostics - Plot of regression diagnostics
%       plotResiduals - Plot of residuals
%       plotSlice - Plot of slices through fitted regression surface
%
%   NonLinearModel properties:
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
%       Robust - Robust fit results
%       Iterative - Summary information about fitting process
%       Formula - Representation of the model used in this fit
%       LogLikelihood - Log of likelihood function at coefficient estimates
%       CoefficientCovariance - Covariance matrix for coefficient estimates
%       CoefficientNames - Coefficient names
%       NumCoefficients - Number of coefficients
%       NumEstimatedCoefficients - Number of estimated coefficients
%
%   See also FITNLM, LinearModel, GeneralizedLinearModel.

%   Copyright 2011-2013 The MathWorks, Inc.

    properties(Constant,Hidden)
        SupportedResidualTypes = {'Raw' 'Pearson' 'Standardized' 'Studentized'};
    end    
    
    properties(GetAccess='public',SetAccess='protected')
%MSE - Mean squared error.
%   The MSE property is the mean squared error. It is an estimate of the
%   variance of the error term.
%
%   See also NonLinearModel.
        MSE = 0;
        
%Iterative - Summary information about fitting process.
%   The Iterative property is a structure with the following fields:
%
%     'InitialCoefs'    Initial coefficient values supplied for the fit
%     'IterOpts'        Options used as input to the fit.
%
%   See also NonLinearModel
        Iterative = [];
        
%Robust - Robust regression results.
%   If the model was fit using robust regression, the Robust property
%   is a structure containing information about that fit. If the model was
%   not fit using robust regression, this property is empty.
%
%   The Robust structure contains the following fields:
%      RobustWgtFun  Weight function used for the fit
%      Tune          Tuning constant used, default [] (selected automatically)
%      Weights       Robust weights used at the final iteration
%
%   See also NonLinearModel.
        Robust = [];
    end
    properties(GetAccess='protected',SetAccess='protected')
        HaveGrad = false;
        Leverage = [];
        Design = []; % not reduced with Subset
    end
    properties(Dependent,GetAccess='public',SetAccess='protected')
%Residuals - Residual values.
%   The Residuals property is a table of residuals. It is a table that
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
%   See also NonLinearModel, plotResiduals, Fitted, predict, random.
        Residuals
                
%Fitted - Fitted (predicted) values.
%   The Fitted property is a vector of fitted values.
%
%   The fitted values are computed using the predictor values used to fit
%   the model. Use the PREDICT method to compute predictions for other
%   predictor values and to compute confidence bounds for the predicted
%   values.
%
%   See also NonLinearModel, Residuals, predict, random.
        Fitted

%RMSE - Root mean squared error.
%   The RMSE property is the root mean squared error. It is an estimate of
%   the standard deviation of the error term.
%
%   See also anova.
        RMSE
        
%Diagnostics - Regression diagnostics.
%   The Diagnostics property is a structure containing a set of diagnostics
%   helpful in finding outliers and influential observations. Many describe
%   the effect on the fit of deleting single observations.  The structure
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
%   CooksDistance is another measure of scaled change in fitted values. A
%   value larger than three times the mean Cook's distance may be
%   considered influential.
%
%   HatMatrix is an N-by-N matrix H such that in a linear model, Yfit=H*Y
%   where Y is the response vector and Yfit is the vector of fitted
%   response values.
%
%   All of these quantities are computed using the usual linear algebra
%   formulas as for a linear regression model, but with the design matrix
%   replaced by the Jacobian of the fitted model.
%
%   See also NonLinearModel.
        Diagnostics
    end
    
    properties(Dependent,GetAccess='public',SetAccess='protected',Hidden=false)    
%WeightedResiduals - Weighted Residual values.
%   The WeightedResiduals property is a table of residuals. It is a table 
%   that has one row for each observation and the following variables:
%
%        'Raw'          Observed minus fitted values
%        'Pearson'      Raw residuals divided by RMSE
%        'Standardized' Raw residuals divided by their estimated standard
%                       deviation
%        'Studentized'  Raw residuals divided by an independent (delete-1)
%                       estimate of their standard deviation
%
%   The ith weighted residual of each type shown above is obtained by multiplying 
%   the corresponding unweighted residual by the square root of effective 
%   observation weight. To obtain any of these columns as a vector, index into 
%   the property using dot notation. For example, in the model M, the weighted 
%   ordinary or weighted raw residual vector is
%
%      r = M.WeightedResiduals.Raw
%
%   See also NonLinearModel, plotResiduals, Fitted, predict, random.
        WeightedResiduals        
    end
    
    properties(Dependent,GetAccess='protected',SetAccess='protected')
        % "Working" values - created and saved during fit, but subject to
        % being cleared out and recreated when needed
        J_r = []; % reduced, i.e. Jac(Subset,:)
    end    
    properties(Access='protected',Hidden=true)        
        %  Structure from nlinfit containing ErrorModel info:
        %
        %  ErrorModel           - Name of the error model.
        %  ErrorParameters      - Parameters of the error model.
        %  ErrorVariance        - Variance function for the error model.
        %  MSE                  - Mean squared error.  
        %  ScheffeSimPred       - Scheffe parameter for simultaneous prediction intervals.
        %  WeightFunction       - True if using 'Weights' with a weight function.
        %  FixedWeights         - True if using 'Weights' with a fixed weight vector.
        %  RobustWeightFunction - True if using Robust fitting with options.RobustWgtFun.
        ErrorModelInfo = struct('ErrorModel',[],'ErrorParameters',[],'ErrorVariance',[],'MSE',[],...
                            'ScheffeSimPred',[],'WeightFunction',[],'FixedWeights',[],'RobustWeightFunction',[]);
        
        % Function handle to a function which returns true if a given contrast
        % of Coefficients is estimable and false otherwise.
        getEstimableContrasts = [];
        
        % Rank of the weighted Jacobian Jw from nlinfit. The coefficient
        % co-variance is mse * inv(Jw'*Jw).
        RankJW = [];
        
        % Store weights supplied as a function handle for use in nlinfit.
        WeightFunctionHandle = [];  
        
        % SST0 - SST0 is the sum of the squared deviations between the observed 
        % response values and 0.
        SST0 = NaN; % total sum of squares, sum(y - 0) over Subset
        
        % F test and p-value for testing the non-linear part of the model.
        % See fTest below for more information.
        FTestFval = NaN;
        FTestPval = NaN;
        EmptyNullModel = NaN;
        HasIntercept = NaN;
        
        % Version number for this class.
        VersionNumber = 2.0;
    end 
    
    properties(Constant=true, Hidden=true)        
         % List of supported error model names.         
            NAME_CONSTANT = 'constant';
        NAME_PROPORTIONAL = 'proportional';
            NAME_COMBINED = 'combined';
     SupportedErrorModels = {NonLinearModel.NAME_CONSTANT,  NonLinearModel.NAME_PROPORTIONAL, NonLinearModel.NAME_COMBINED};          
    end
    
    methods % get/set       
        
        function yfit = get.Fitted(model)
            yfit = get_fitted(model);
        end
        function r = get.Residuals(model)
            r = get_residuals(model);            
        end
        function r = get.WeightedResiduals(model)
            r = get_weighted_residuals(model);            
        end
         function s = get.RMSE(model)
            s = sqrt(model.MSE);
        end
        function J_r = get.J_r(model)
            if isempty(model.WorkingValues)
                J_r = create_J_r(model);
            else
                J_r = model.WorkingValues.J_r;
            end
        end
        
% The following code is removed because it causes a bad interaction with
% the array editor. As a result, the Diagnostics propety does not appear in
% the array editor view of the LinearModel object. Diagnostics property
% access from the command line is provided in the subsref method.

%         function D = get.Diagnostics(model)
%             D = get_diagnostics(model);
%         end
    end % get/set
    
    methods(Hidden=true, Access='public')
        function model = NonLinearModel(varargin) % modeldef,coefs,etc.
            if nargin == 0 % special case
                model.Formula = classreg.regr.NonLinearFormula;
                return
            end
            error(message('stats:NonLinearModel:NoConstructor'));
        end
    end
    
    methods(Access='public')
        function disp(model)
%DISP Display a NonLinearModel.
%   DISP(NLM) displays the NonLinearModel NLM.
%
%   See also NonLinearModel.
            isLoose = strcmp(get(0,'FormatSpacing'),'loose');
            if (isLoose), fprintf('\n'); end
            if isempty(model.Robust)
                fprintf(getString(message('stats:NonLinearModel:display_Nonlin')));
            else
                fprintf(getString(message('stats:NonLinearModel:display_NonlinRobust')));
            end

            % Fit the formula string to the command window width
            indent = '    ';
            maxWidth = matlab.desktop.commandwindow.size; maxWidth = maxWidth(1) - 1;
            f = model.Formula;
            fstr = char(f,maxWidth-length(indent));
            disp([indent fstr]);
            
            if model.IsFitFromData
                fprintf(getString(message('stats:NonLinearModel:display_EstimatedCoefficients')));
                disp(model.Coefficients);
                fprintf('%s',getString(message('stats:NonLinearModel:display_NumObsDFE',model.NumObservations,model.DFE)));
                fprintf('%s',getString(message('stats:NonLinearModel:display_RMSE',num2str(model.RMSE,'%.3g'))));
                rsq = get_rsquared(model,{'ordinary','adjusted'});
                fprintf('%s',getString(message('stats:NonLinearModel:display_RSquared',num2str(rsq(1),'%.3g'),num2str(rsq(2),'%.3g'))));
                [f,p, emptyNullModel] = fTest(model);
                if emptyNullModel
                    fprintf('%s',getString(message('stats:NonLinearModel:display_FtestZero',num2str(f,'%.3g'),num2str(p,'%.3g'))));
                else
                    fprintf('%s',getString(message('stats:NonLinearModel:display_Ftest',num2str(f,'%.3g'),num2str(p,'%.3g'))));
                end
            else
                fprintf(getString(message('stats:NonLinearModel:display_Coefficients')));
                if any(model.Coefficients.SE > 0)
                    disp(model.Coefficients(:,{'Value' 'SE'}));
                else
                    disp(model.Coefficients(:,{'Value'}));
                end
                if model.MSE > 0
                    fprintf('\n%s',getString(message('stats:NonLinearModel:display_RMSE',num2str(model.RMSE,'%.3g'))));
                end
            end
        end
        
        % --------------------------------------------------------------------
        function [ypred,yCI] = predict(model,varargin)
%PREDICT Compute predicted values given predictor values.
%   YPRED = PREDICT(NLM,DS) computes a vector YPRED of predicted values from
%   the NonLinearModel NLM using predictor variables from the dataset/table DS. DS
%   must contain all of the predictor variables used to create NLM.
%
%   YPRED = PREDICT(NLM,X), where X is a data matrix with the same number of
%   columns as the matrix used to create NLM, computes predictions using the
%   values in X.
%
%   [YPRED,YCI] = PREDICT(...) also returns the two-column matrix YCI
%   containing 95% confidence intervals for the predicted values. These are
%   non-simultaneous intervals for predicting the mean response at the
%   specified predictor values. The lower limits of the bounds are in
%   column 1, and the upper limits are in column 2.
%
%   [...] = PREDICT(NLM,DS,PARAM1,VAL1,PARAM2,VAL2,...) or
%   [...] = PREDICT(NLM,X,PARAM1,VAL1,PARAM2,VAL2,...) specifies one or more
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
%      'Weights'      A vector W of real positive values with the same number
%                     of elements as the number of observations (or rows) in 
%                     X. The error variance at observation i is estimated as 
%                     MSE * (1/W(i)) where MSE is the mean squared error.               
%                     'Weights' can also be specified as a function handle that 
%                     accepts a vector of predicted response values and returns 
%                     a vector of real positive values as output. Default is no 
%                     weights.
%
%   Example:
%      % Compute predicted values and confidence bounds for three
%      % observations
%      load reaction
%      nlm = fitnlm(reactants,rate,@hougen,[1 .05 .02 .1 2])
%      [yfit,confint] = predict(nlm,reactants(1:3,:))
%
%   See also NonLinearModel, random.

            if nargin > 1 && ~internal.stats.isString(varargin{1})
                Xpred = varargin{1};
                varargin = varargin(2:end);
                design = predictorMatrix(model,Xpred);
            else
                design = model.Design;
            end
            
            paramNames = {'Confidence', 'Simultaneous', 'Prediction', 'Alpha', 'Weights'};
            paramDflts = {        .95 ,         false ,      'curve',  .05   ,     []   };
            [conf,simOpt,predOpt,alpha,weights,supplied] = ...
                internal.stats.parseArgs(paramNames, paramDflts, varargin{:});
            
            if supplied.Confidence && supplied.Alpha
                error(message('stats:NonLinearModel:ArgCombination'))
            end
            if supplied.Alpha
                conf = 1-alpha;
            end
            if ~islogical(predOpt)
                predOpt = internal.stats.getParamVal(predOpt,{'curve' 'observation'},'''Prediction''');
                predOpt = isequal(predOpt,'observation');
            end
            if nargout < 2
                ypred = model.Formula.ModelFun(model.Coefs,design);
            else
                % Recomputing the Jacobian here will allow CIs for excluded
                % observations, and give NaN CIs for missing observations
                if model.HaveGrad
                    [ypred,Jpred] = model.Formula.ModelFun(model.Coefs,design);
                else
                    ypred = model.Formula.ModelFun(model.Coefs,design);
                    Jpred = jacobian(model,design);
                end
                
                % (1) Which rows of Jpred are estimable?
                estimable = model.getEstimableContrasts(Jpred);
                
                % (2) Make sure weights are sensible.
                if ~isempty(weights)
                    weights = NonLinearModel.checkWeights(weights, ypred);
                end
                
                % (3) Compute variance of the predictions.
                ypredVar = sum((Jpred*model.CoefficientCovariance) .* Jpred,2);
                % Are we using an error model to compute errorVar?
                usingErrorModel = false;
                if (predOpt)
                    if ~isempty(weights)
                        if isa(weights,'function_handle')
                            wVec = weights(ypred);
                        else
                            wVec = weights;
                        end
                        % Assume error variance is mse*(1./wVec).
                        wVec = wVec(:);
                        errorVar = model.MSE./wVec;
                    elseif ~isempty(model.ErrorModelInfo) && ~isempty(model.ErrorModelInfo.ScheffeSimPred) && ~isempty(model.ErrorModelInfo.ErrorVariance)
                        % Use ErrorVariance function in model.ErrorModelInfo.
                        errorVar = model.ErrorModelInfo.ErrorVariance(design);
                        errorVar = errorVar(:);
                        usingErrorModel = true;
                    else
                        % Assume a constant variance model if errorModelInfo and weights are 
                        % not supplied.
                        errorVar = model.MSE * ones(size(Jpred,1),1);
                    end    
                    ypredVar = ypredVar + errorVar;
                end                
                                              
                % (4) Compute delta, the CI half-widths.
                if (simOpt)
                   % Simultaneous CIs.
                   if (predOpt)
                        % Simultaneous CI on new observations.
                        if usingErrorModel
                            % Use the precomputed Scheffe parameter valid for using this error model.
                            sch = model.ErrorModelInfo.ScheffeSimPred;
                            % Make sure sch is either rankJ or (rankJ+1).
                            if ( sch ~= model.RankJW && sch ~= (model.RankJW+1) )
                                % Be conservative.
                                sch = (model.RankJW+1);
                            end            
                            % Validate Scheffe parameter against Jacobian at prediction points.
                            sch = NonLinearModel.validateScheffeParamUsingJacobianPred(sch, model.RankJW, Jpred, estimable, errorVar);                                    
                        else
                            % Be conservative.
                            sch = (model.RankJW+1);
                        end        
                   else
                       % Simultaneous CI on curve.
                        sch = model.RankJW;
                   end
                   crit = sqrt(sch * finv(conf, sch, model.DFE));   
                else
                   % Pointwise CIs.
                   crit = tinv((1+conf)/2,model.DFE);
                end
                delta = sqrt(ypredVar) * crit;
                
                % (5) Compute CIs.
                yCI = [ypred-delta, ypred+delta];               
                % Set CIs for non-estimable rows of Jpred to NaN.
                yCI(~estimable,:) = NaN;
            end
        end
        
        % --------------------------------------------------------------------
        function ysim = random(model,varargin) % model, Xsim
%RANDOM Generate random response values given predictor values.
%   YNEW = RANDOM(NLM,DS) generates a vector YNEW of random values from the
%   NonLinearModel NLM using predictor variables from the dataset/table DS. DS
%   must contain all of the predictor variables used to create NLM. The
%   output YNEW is computed by creating predicted values and adding new
%   random noise values with standard deviation NLM.RMSE.
%
%   YNEW = RANDOM(NLM,X), where X is a data matrix with the same number of
%   columns as the matrix used to create NLM, generates responses using the
%   values in X.
%
%   [...] = RANDOM(NLM,DS,PARAM1,VAL1,PARAM2,VAL2,...) or
%   [...] = RANDOM(NLM,X,PARAM1,VAL1,PARAM2,VAL2,...) allows you to pass an 
%   optional name/value pair that specifies the observation weights: 
%
%      'Weights'      A vector W of real positive values with the same number
%                     of elements as the number of observations (or rows) in 
%                     X. The error variance at observation i is estimated as 
%                     MSE * (1/W(i)) where MSE is the mean squared error.               
%                     'Weights' can also be specified as a function handle that 
%                     accepts a vector of predicted response values and returns 
%                     a vector of real positive values as output. Default is no 
%                     weights.
%
%   Example:
%      % Fit a nonlinear model to census data
%      load census
%      start = [1 1 1900 1];
%      nlm = fitnlm(cdate,pop,'pop ~ b1 + b2*abs(cdate-b3)^b4',start)
%
%      % See if the real data look like a random sample from this model
%      subplot(1,2,1)
%      plot(cdate,pop,'bo',cdate,predict(nlm,cdate),'r-')
%
%      pop2 = random(nlm,cdate);
%      subplot(1,2,2)
%      plot(cdate,pop2,'bo',cdate,predict(nlm,cdate),'r-')
%
%   See also NonLinearModel, predict, RMSE.


            if nargin > 1 && ~internal.stats.isString(varargin{1})
                Xpred = varargin{1};
                varargin = varargin(2:end);
                design = predictorMatrix(model,Xpred);
            else
                design = model.Design;
            end
            
            paramNames = {'Weights'};
            paramDflts = {   []    };
            [weights,~] = internal.stats.parseArgs(paramNames, paramDflts, varargin{:});
            
            % Get ypred.
            ypred = model.Formula.ModelFun(model.Coefs,design);

            % (1) Make sure weights are sensible.
            if ~isempty(weights)
                weights = NonLinearModel.checkWeights(weights, ypred);
            end
            % (2) Compute variance of the predictions.
            if ~isempty(weights)
                if isa(weights,'function_handle')
                    wVec = weights(ypred);
                else
                    wVec = weights;
                end
                % Assume error variance is mse*(1./wVec).
                wVec = wVec(:);
                errorVar = model.MSE./wVec;
            elseif ~isempty(model.ErrorModelInfo) && ~isempty(model.ErrorModelInfo.ScheffeSimPred) && ~isempty(model.ErrorModelInfo.ErrorVariance)
                % Use ErrorVariance function in model.ErrorModelInfo.
                errorVar = model.ErrorModelInfo.ErrorVariance(design);
                errorVar = errorVar(:);
            else
                % Assume a constant variance model if errorModelInfo and weights are 
                % not supplied.
                errorVar = model.MSE * ones(length(ypred),1);
            end    
            
            % (3) Assume error distribution is Gaussian for simulation
            % purposes. We did not make this assumption when using GLS
            % for model fitting.
            ysim = normrnd(ypred,sqrt(errorVar));
        end
        
        % -------------------- pass throughs to modelutils -------------------
        function hout = plotDiagnostics(model,plottype,varargin)
%plotDiagnostics Plot diagnostics of fitted model
%    plotDiagnostics(NLM,PLOTTYPE) plots diagnostics from NonLinearModel NLM
%    in a plot of type PLOTTYPE. The default value of PLOTTYPE is 'leverage'.
%    Valid values for PLOTTYPE are:
%
%       'contour'      residual vs. leverage with overlaid Cook's contours
%       'cookd'        Cook's distance
%       'leverage'     leverage (diagonal of Hat matrix)
%
%    H = plotDiagnostics(...) returns handles to the lines in the plot.
%
%    The PLOTTYPE argument can be followed by parameter/value pairs to
%    specify additional properties of the primary line in the plot. For
%    example, plotDiagnostics(NLM,'cookd','Marker','s') uses a square
%    marker.
%
%    The data cursor tool in the figure window will display the X and Y
%    values for any data point, along with the observation name or number.
%    It also displays the coefficient name for 'dfbetas'.
%
%   Example:
%      % Fit the Hougen model to reaction data, and identify points with
%      % high leverage
%      load reaction
%      nlm = fitnlm(reactants,rate,@hougen,[1 .05 .02 .1 2])
%      plotDiagnostics(nlm,'leverage')
%      find(nlm.Diagnostics.Leverage > 0.8)
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
%        'Raw'          Observed minus fitted values
%        'Pearson'      Raw residuals divided by RMSE
%        'Standardized' Raw residuals divided by their estimated standard
%                       deviation
%        'Studentized'  Raw residuals divided by an independent (delete-1)
%                       estimate of their standard deviation
%
%    For the Standardized and Studentized choices, the calculations are
%    based on a linear approximation to the nonlinear function at the
%    estimated coefficient values.
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
%   Example:
%      % Fit the Hougen model to reaction data, and make a normal
%      % probability plot of the residuals
%      load reaction
%      nlm = fitnlm(reactants,rate,@hougen,[1 .05 .02 .1 2])
%      plotResiduals(nlm,'probability')
%
%    See also NonLinearModel, plotDiagnostics.
            if nargin<2
                plottype = 'histogram';
            end
            [residtype,wantWeighted,~,args] = internal.stats.parseArgs({'residualtype','weighted'},{'Raw',false},varargin{:});
            varargin = args;
            residtype = internal.stats.getParamVal(residtype,...
                         NonLinearModel.SupportedResidualTypes,'''ResidualType''');            
            wantWeighted = internal.stats.parseOnOff(wantWeighted,'''Weighted''');
            f = classreg.regr.modelutils.plotResiduals(model,plottype,'ResidualType',residtype,'Weighted',wantWeighted,varargin{:});
            if nargout>0
                hout = f;
            end
        end
        
        function fout = plotSlice(model)
%plotSlice Plot slices through the fitted regression surface.
%   plotSlice(NLM) creates a new figure containing a series of plots, each
%   representing a slice through the regression surface predicted by NLM.
%   For each plot, the surface slice is shown as a function of a single
%   predictor variable, with the other predictor variables held constant.
%
%   Example:
%      % Fit the Hougen model to reaction data, and make a slice
%      % plot of the fitted regression surface
%      load reaction
%      vnames = {'Hydrogen' 'nPentane' 'Isopentane' 'ReactionRate'};
%      nlm = fitnlm(reactants,rate,@hougen,...
%                               [1 .05 .02 .1 2],'VarNames',vnames)
%      plotSlice(nlm)
%
%   See also NonLinearModel, predict, RMSE.
            f = classreg.regr.modelutils.plotSlice(model);
            if nargout>0
                fout = f;
            end
        end
    end % public methods
    
    methods(Access='public')
        function CI = coefCI(model,alpha) % a method, to accept alpha argument
 %coefCI Confidence intervals for coefficients.
%   CI = coefCI(M) computes 95% confidence intervals for the coefficients
%   in the regression model M. The output CI is a two-column matrix with
%   the lower confidence limits in column 1 and the upper confidence limits
%   in column 2.
%
%   CI = coefCI(M,ALPHA) computes 100*(1-ALPHA)% confidence intervals. The
%   default is ALPHA=0.05 for 95% confidence.
%
%   Example:
%       % Find 90% confidence limits for estimated coefficients
%       load reaction
%       nlm = fitnlm(reactants,rate,@hougen,[1 .05 .02 .1 2])
%       confint = coefCI(nlm,0.1)
%
%   See also LinearModel, GeneralizedLinearModel, NonLinearModel, linhyptest.
           if nargin < 2
                alpha = 0.05;
            end
            se = sqrt(diag(model.CoefficientCovariance));
            delta = se * tinv(1-alpha/2,model.DFE);
            CI = [(model.Coefs(:) - delta) (model.Coefs(:) + delta)];
            
            % Make CIs corresponding to non-estimable coefficients equal to [NaN, NaN].
            nC = numel(model.Coefs);
            H = eye(nC);
            estimable = model.getEstimableContrasts(H);
            CI(~estimable,:) = NaN;
        end
        
        % --------------------------------------------------------------------
        function [p,t,r] = coefTest(model,H,c)
%coefTest Linear hypothesis test on coefficients.
%   P = coefTest(M) computes the p-value for an F test that all
%   coefficient estimates in the regression model M are zero.
%
%   P = coefTest(M,H), with H a numeric matrix having one column for each
%   coefficient, performs an F test that H*B=0, where B represents the
%   coefficient vector.
%
%   P = coefTest(M,H,C) accepts a vector C having one value for each row
%   of H, and it performs an F test that H*B=C.
%
%   [P,F,R] = coefTest(...) also returns the F-statistic F and the rank R
%   of the matrix H. The F statistic has R degrees of freedom in the
%   numerator and M.DFE degrees of freedom in the denominator.
%
%   Example:
%      % Test to see if all coefficients in the denominator of the
%      % model expression could be the same
%      load reaction
%      myfun = 'rate~(b1*x2-x3/b5)/(1+b2*x1+b3*x2+b4*x3)';
%      nlm = fitnlm(reactants,rate,myfun,[1 .05 .02 .1 2])
%      p = coefTest(nlm,[0 1 -1 0 0;0 0 1 -1 0])
%
%   See also LinearModel, GeneralizedLinearModel, NonLinearModel, linhyptest.

            % Test that Hb == c where b is the coefficient vector.
            nc = model.NumCoefficients;
            if nargin<2
                % Default is to test all terms = 0
                H = eye(nc);
            end
            if nargin < 3
                c = zeros(size(H,1),1);
            end
            % Outputs p = p-value, t = F statistic, r = rank of test.
            [p,t,r] = linhyptest(model.Coefs,model.CoefficientCovariance,c,H,model.DFE);
            
            % If all rows of H are not estimable, return NaNs.
            estimable = model.getEstimableContrasts(H);
            if ~all(estimable)
                p = p*NaN;
                t = t*NaN;
                r = r*NaN;
            end
        end
    end % public methods
    
    methods(Hidden=true)
        
        function t = title(model)
            % Fit the formula string to the command window width
            indent = '    ';
            maxWidth = matlab.desktop.commandwindow.size; maxWidth = maxWidth(1) - 1;
            f = model.Formula;
            fstr = char(f,maxWidth-length(indent));
            t = [indent fstr];
        end
        
        function v = varianceParam(model)
            v = model.RMSE;
        end
    end        

    methods(Access='protected')
                          
        % Test if non-linear part of the model function is zero. Do one of
        % the following:
        % (1) If full model has an intercept and number of estimated coefficients
        %     is > 1 then test full model against the intercept only model.  
        % (2) Otherwise, test full model against a zero model.
        function [f,p,emptyNullModel,hasIntercept] = fTest(model)
        % Variable emptyNullModel is:
        %  -  true if our null model is the zero model. 
        %  - false if our null model is the intercept only model.            
            
            % (1) Inspect unweighted Jacobian and figure out if there is an
            % intercept (i.e., if if Junw_r has a constant column).
            [~,Junw_r] = create_J_r(model);
            Jmin = min(Junw_r,[],1);
            Jmax = max(Junw_r,[],1);
            hasIntercept = any( abs(Jmax-Jmin) <= sqrt(eps(class(Junw_r))) * (abs(Jmax) + abs(Jmin)) ); 
            if hasIntercept &&  (model.NumEstimatedCoefficients > 1)
                % (2) Compare full model vs. intercept only model.
                emptyNullModel = false;
                nobs = model.NumObservations;
                ssr = max(model.SST - model.SSE,0);
                dfr = model.NumEstimatedCoefficients - 1;
                dfe = nobs - 1 - dfr;
                f = (ssr./dfr) / (model.SSE/dfe);
                p = fcdf(1./f,dfe,dfr); % upper tail
            else
                % (2) Compare full model vs. zero model.
            	emptyNullModel = true;
                ssr = max(model.SST0 - model.SSE,0);
                dfr = model.NumEstimatedCoefficients;
                dfe = model.NumObservations - model.NumEstimatedCoefficients;
                f = (ssr./dfr) / (model.SSE/dfe);
                p = fcdf(1./f,dfe,dfr); % upper tail
            end                        
        end        
        
        % Modify the standard display table to display only estimable
        % coefficients.
        function tbl = tstats(model)
            
            % Account for rank deficient case before calling tstats.
            effectiveNumObservations = model.DFE + numel(model.CoefficientNames);            
            tbl = classreg.regr.modelutils.tstats(model.Coefs,sqrt(diag(model.CoefficientCovariance)), ...
                                      effectiveNumObservations,model.CoefficientNames);
                                                                                                     
            % (1) Get coefficients that are estimable.
            % (2) Get the contents of table and replace columns 2-4 
            %     of rows for non-estimable coefficients by NaNs.
            % (3) Put the NaN filled tableVals back into the original table.
            estimable = model.getEstimableContrasts(eye(numel(model.Coefs)));            
            tableVals = table2array(tbl);
            tableVals(~estimable,2:4) = NaN;
            tableVals = array2table(tableVals);
            tableVals.Properties.VariableNames = tbl.Properties.VariableNames;
            tableVals.Properties.RowNames = tbl.Properties.RowNames;
        end
        
       % "overloaded" get functions
        function D = get_diagnostics(model,type)
            if nargin<2 % return all diagnostics in a table
                HatMatrix = get_diagnostics(model,'hatmatrix');
                CooksDistance = get_diagnostics(model,'cooksdistance');
                Leverage = model.Leverage;
                D = table(Leverage,CooksDistance,HatMatrix,...
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
                        % D should already satisfy this: D(~subset,~subset) = 0;                  
                    case 'cooksdistance'
                        D = get_CooksDistance(model);
                        D(~subset,:) = NaN;
                    otherwise
                        error(message('stats:LinearModel:UnrecognizedDiagnostic', type));
                end
            end
        end
        function s2_i = get_S2_i(model)
            r = getResponse(model) - predict(model);
            h = model.Leverage;
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
            subset = model.ObservationInfo.Subset;
            SSE = sum(wt(subset) .* abs(r(subset)).^2);
            s2_i = (SSE - delta_i) ./ newdf;
            s2_i(~subset & ~isnan(s2_i)) = 0;
        end

        function H = get_HatMatrix(model)
            if hasData(model)
                % Let Xw = sqrt(W)*X.
                % yfit = X * inv(X'*W*X) * X'*W*y 
                %      = sqrt(W^-1)* [Xw*inv(Xw'*Xw)*Xw] *sqrt(W)*y
                %      = H*y
                % Compute [...] using a QR decomposition on Xw.
                
                % Get the unweighted Jacobian, because the weighted version
                % has only the prior weights and we also need robust
                % weights
                [~,Xunw_r] = create_J_r(model);

                sw = sqrt(get_CombinedWeights_r(model));
                Xw_r = bsxfun(@times, sw(:), Xunw_r);
                [Qw,~,~] = qr(Xw_r,0); % do the same pivoting as in the real fit
                rank = model.NumEstimatedCoefficients;
                Qw = Qw(:,1:rank);
                T = Qw*Qw';
                H = zeros(model.NumObservations);
                subset = model.ObservationInfo.Subset;
                H1 = bsxfun(@times,bsxfun(@times,1./sw,T),sw');
                H1(sw<=0,:) = 0; % prefer 0 to NaN for pt with 0 weight
                H(subset,subset) = H1;
            else
                H = [];
            end
        end
        
        function d = get_CooksDistance(model)
            if hasData(model)
                w = get_CombinedWeights_r(model,false);
                r = model.Residuals.Raw;
                h = model.Leverage;
                d = w .* abs(r).^2 .* (h./((1-h).^2))./(model.NumCoefficients*model.MSE);
            else
                d = [];
            end
        end
        
        function w = get_CombinedWeights_r(model,reduce)
            % Combined weights are obtained as follows:
            %
            % (1) For a fixed weight vector (possibly all ones) and NO
            %     robust regression, just return model.ObservationInfo.Weights.
            %     Make sure that model.Robust is [] using an assertion.
            %
            % (2) For a fixed weight vector (possibly all ones) and
            %     using robust regression, return the product of
            %     model.ObservationInfo.Weights and model.Robust.Weights.
            %     Since we currently disallow weight vector and robust options 
            %     simultaneously, model.Observation.Weights must be all ones.
            %     Hence the assertion below.
            %
            % (3) When using weights as a function handle or when using the
            %     combined or proportional error model, return an effective
            %     weight given by:
            %       w = model.MSE ./ model.ErrorModelInfo.ErrorVariance(X) 
            %     where X = model.Design.
            %
            %  For reference, here are the fields in ErrorModelInfo:
            %
            %  ErrorModel           - Name of the error model.
            %  ErrorParameters      - Parameters of the error model.
            %  ErrorVariance        - Variance function for the error model.
            %  MSE                  - Mean squared error.  
            %  ScheffeSimPred       - Scheffe parameter for simultaneous prediction intervals.
            %  WeightFunction       - True if using 'Weights' with a weight function.
            %  FixedWeights         - True if using 'Weights' with a fixed weight vector.
            %  RobustWeightFunction - True if using Robust fitting with options.RobustWgtFun.
            
            switch lower(model.ErrorModelInfo.ErrorModel)
                case NonLinearModel.NAME_CONSTANT                    
                    if ~isempty(model.WeightFunctionHandle)
                        % Weights as a function handle. Get effective weight.
                        X = model.Design;
                        w = model.MSE ./ model.ErrorModelInfo.ErrorVariance(X);
                        % Assert that model.Robust is [].
                        assert(isempty(model.Robust));
                    elseif ~isempty(model.Robust)
                        % Robust weight function. 
                        w = model.ObservationInfo.Weights .* model.Robust.Weights;
                        % Since we disallow weight vector and robust options 
                        % simultaneously, model.Observation.Weights must be all ones.
                        assert(all(model.ObservationInfo.Weights==1));
                    else 
                        % Fixed weight vector or a weight vector of all ones.
                        w = model.ObservationInfo.Weights;
                        % Assert that model.Robust is [].
                        assert(isempty(model.Robust));
                    end                                                           
                case {NonLinearModel.NAME_COMBINED,NonLinearModel.NAME_PROPORTIONAL}
                    % Non constant error model. Get effective weight.
                    X = model.Design;
                    w = model.MSE ./ model.ErrorModelInfo.ErrorVariance(X);                
                otherwise
                    error(message('stats:NonLinearModel:InvalidErrorModel'));
            end
    
            if nargin<2 || reduce
                subset = model.ObservationInfo.Subset;
                w = w(subset);
            end
        end
       
        % -------------------- pass throughs to modelutils -------------------
        function design = predictorMatrix(model,X,isPredictorsOnly)
            if isa(X,'dataset')
                X = dataset2table(X);
            end
            if isa(X,'table')              
                % Given a dataset/table, make sure it has all of the predictor
                % variables.  It can also have any unused variables from
                % the model fit.
                tf = ismember(model.PredictorNames,X.Properties.VariableNames);
                if ~all(tf)
                    error(message('stats:NonLinearModel:PredictorMissing'));
                end
                design = classreg.regr.modelutils.predictormatrix(X,...
                    'ResponseVar',[],...
                    'PredictorVars',model.Formula.PredictorNames);
            else
                if nargin > 2 && isPredictorsOnly
                    % Given a matrix, assume the columns correspond in order
                    % to the predictor variables
                    npreds = model.NumPredictors;
                    if size(X,2) ~= npreds
                        error(message('stats:NonLinearModel:BadXColumnNumber', npreds));
                    end
                    varLocs = model.PredLocs;
                else
                    % Given a matrix, assume the columns correspond in order
                    % to the non-response variables
                    nvars = model.NumVariables;
                    if size(X,2) ~= nvars-1
                        error(message('stats:NonLinearModel:BadXColumnNumber', nvars - 1));
                    end
                    varLocs = true(1,nvars);
                    varLocs(model.RespLoc) = false;
                end
                design = X(:,model.VariableInfo.InModel(varLocs));
            end
        end
            
        % --------------------------------------------------------------------
        function J = jacobian(model,X)
            if model.HaveGrad
                [~,J] = model.Formula.ModelFun(model.Coefs, X);
            else
                % Approximate the Jacobian at Formula.ModelFun(beta,X).
                beta = model.Coefs;
                dbeta = eps(max(abs(beta),1)).^(1/3);
                J = zeros(size(X,1),model.NumCoefficients);
                for i = 1:model.NumCoefficients
                   h = zeros(size(beta)); h(i) = dbeta(i);
                   ypredplus = model.Formula.ModelFun(beta+h, X);
                   ypredminus = model.Formula.ModelFun(beta-h, X);
                   J(:,i) = (ypredplus - ypredminus) / (2*h(i));
                end
            end
        end
        
        % --------------------------------------------------------------------
        function model = fitter(model)
            
            % Get Full predictor matrix and Full response.
            X = getData(model);
            model.Design = predictorMatrix(model,X);
            model.CoefficientNames = model.Formula.CoefficientNames;
            response = getResponse(model);

            % Get everything needed to call nlinfit and use selected
            % observations as per model.ObservationInfo.Subset.
            opts = model.Iterative.IterOpts;
            subset = model.ObservationInfo.Subset;
            X = model.Design(subset,:);
            y = response(subset,:);
            F = model.Formula.ModelFun;
            b0 = model.Iterative.InitialCoefs;
            
            
            % For weighted fits we need to minimize sum(w.*(y-f).^2). But
            % we can absorb sqrt(w) into both y and f. So define a new
            % response vector and a new model function. All this is done
            % inside nlinfit. Create param/value pair for argument
            % 'Weights' which may have been specified as a function handle.
            if isempty(model.WeightFunctionHandle)
                % Vector of weights.
                w = model.ObservationInfo.Weights(subset,:);
                if ~all(w==1)
                    wtargs = {'Weights',w};
                else
                    wtargs = {};
                end
            else
                % Weights = function handle.
                w = model.WeightFunctionHandle;
                wtargs = {'Weights',w};
            end
            
            % Create param/value pair for arguments 'ErrorModel' and 'ErrorModelParameters'.
            errormodelargs = {'ErrorModel',model.ErrorModelInfo.ErrorModel,...
                    'ErrorParam',model.ErrorModelInfo.ErrorParam};
            
            % Both Weights and ErrorModels cannot be used unless ErrorModel 
            % is 'constant'.
            if ~isempty(wtargs) && ~strcmpi(model.ErrorModelInfo.ErrorModel,NonLinearModel.NAME_CONSTANT)
                error(message('stats:NonLinearModel:ErrorModelWeightConflict'));
            end
            
            % Call nlinfit and collect outputs. Here's the sequence of outputs 
            % from nlinfit: [beta,r,J,Sigma,mse,errorModelInfo,robustw].
            if isempty(model.Robust)
                [model.Coefs,~,J_r,model.CoefficientCovariance,model.MSE,model.ErrorModelInfo,~] = ...
                    nlinfit(X,y,F,b0,opts,wtargs{:},errormodelargs{:});
            else
                opts.Robust = 'on';
                opts = statset(opts,model.Robust);
                [model.Coefs,~,J_r,model.CoefficientCovariance,model.MSE,model.ErrorModelInfo,robustw] = ...
                    nlinfit(X,y,F,b0,opts,wtargs{:},errormodelargs{:});
                % Set robust weights from final iteration of robust fitting 
                % but for those observations that were actually used in fitting.
                model.Robust.Weights = zeros(length(subset),1);
                model.Robust.Weights(subset) = robustw;
            end
            
            % Populate the J_r field in the WorkingValues structure.
            model.WorkingValues.J_r = J_r;
            
            % Create a function handle to return estimable contrasts (used 
            % in predict and fTest) and compute the rank of J_r (used to compute model.DFE). 
            TolSVD = eps(class(model.Coefs));
            [~,model.RankJW,~,model.getEstimableContrasts,~,~] = ...
                internal.stats.isEstimable(eye(numel(model.Coefs)),'DesignMatrix',J_r,'TolSVD',TolSVD);
            
            % Get effective degrees of freedom.
            model.DFE = model.NumObservations - model.RankJW;                                  
        end
            
        % --------------------------------------------------------------------
        function model = selectVariables(model)
            f = model.Formula;
            [~,model.PredLocs] = ismember(f.PredictorNames,f.VariableNames);
            [~,model.RespLoc] = ismember(f.ResponseName,f.VariableNames);
            model = selectVariables@classreg.regr.ParametricRegression(model);
        end
        
        % --------------------------------------------------------------------
        function model = postFit(model)           
            %model = postFit@classreg.regr.ParametricRegression(model);
            model = getSumOfSquaresAndLogLikelihood(model);
            [~,Junw_r] = create_J_r(model);                        
            sw = sqrt(get_CombinedWeights_r(model));
            Xw_r = bsxfun(@times,Junw_r,sw);
            [Qw_r,~,~] = qr(Xw_r,0); % do the same pivoting as in the real fit
            
            % model.NumEstimatedCoefficients = model.NumObservations - model.DFE.
            % model.NumObservations is the observations actually used in the fit.
            % Since model.DFE is correctly set in method fitter, rank should be 
            % correct. Also note that QR factorization with 3 output args orders 
            % the diagonal elements of R in decreasing absolute value. Hence, we 
            % can just use the 1:rank columns of Qw_r.
            rank = model.NumEstimatedCoefficients;
            Qw_r = Qw_r(:,1:rank);
            h = zeros(size(model.ObservationInfo,1),1);
            h(model.ObservationInfo.Subset) = sum(abs(Qw_r).^2,2);
            model.Leverage = h;
            
            % Test if non-linear part of the model function is zero.
            [model.FTestFval,model.FTestPval,model.EmptyNullModel,...
                model.HasIntercept] = fTest(model);           
        end
        
        %---------------------------------------------------------------------
        function model = getSumOfSquaresAndLogLikelihood(model)
            % Convenience function to get weighted sum of squares and log
            % likelihoods. Note that we use effective weights returned by
            % get_CombinedWeights_r - not observation weights.
            
            % Subset of observations.
            subset = model.ObservationInfo.Subset;
            % Subset of: raw residuals = resid_r = get_residuals(model,'raw');
            resid_r = getResponse(model) - predict(model);
            resid_r = resid_r(subset);
            % Subset of response.
            y = getResponse(model);
            y_r = y(subset);
            % Subset of model fit.
            yfit_r = predict(model);
            yfit_r = yfit_r(subset);
            % Subset of combined weights.
            w = get_CombinedWeights_r(model,false);
            w_r = w(subset);
            % Calculations for Sum of Squares and Log Likelihoods.
            sumw = sum(w_r);
            wtdymean = sum(w_r.*y_r) / sumw;
            model.SSE = sum(w_r .* resid_r.^2);
            model.SSR = sum(w_r .* (yfit_r - wtdymean).^2);
            model.SST = sum(w_r .* (y_r - wtdymean).^2);
           model.SST0 = sum(w_r .* (y_r).^2);
            model.LogLikelihood = getlogLikelihood(model);
            model.LogLikelihoodNull = logLikelihoodNull(model);            
        end
        
        % --------------------------------------------------------------------
        function [J_r,Junw_r] = create_J_r(model)
            subset = model.ObservationInfo.Subset;
            J_r = jacobian(model,model.Design(subset,:));
            if nargout>=2
                Junw_r = J_r; % unweighted version of Jacobian
            end
            w = model.ObservationInfo.Weights(subset,:);
            if ~all(w==1)
                J_r = bsxfun(@times, sqrt(w), J_r);
            end
        end
        
        % --------------------------------------------------------------------
        function ypred = predictPredictorMatrix(model,Xpred)
            % Assume the columns correspond in order to the predictor variables
            design = predictorMatrix(model,Xpred,true);
            ypred = model.Formula.ModelFun(model.Coefs,design);
        end
        
        % --------------------------------------------------------------------
        function L = getlogLikelihood(model)
            Var = model.DFE/model.NumObservations*model.MSE;
            if isempty(model.Robust)
                % Get effective weights for observations used in the fit.            
                w_r = get_CombinedWeights_r(model,true);
                L = -(model.DFE + model.NumObservations*log(2*pi) + sum(log(Var./w_r)))/2;
            else
                % Return the Gaussian log likelihood using the robust
                % estimate of MSE.
                subset = model.ObservationInfo.Subset;
                yfit = predict(model);
                yfit = yfit(subset);
                y = getResponse(model);
                y = y(subset);
                L = -( sum((y-yfit).^2)/Var + model.NumObservations*log(2*pi) + model.NumObservations*log(Var))/2;
            end
            
        end
        
        % --------------------------------------------------------------------
        function L0 = logLikelihoodNull(model)
            
            % A null model is a constant only model. So the model function
            % looks like: linfun0 = @(beta,X) beta * ones(size(X,1),1) 
            % and beta is a scalar parameter. The error variance is estimated 
            % as model.MSE for Robust fitting and model.MSE/w(i) for
            % observation i when using non-Robust fitting. For non-Robust 
            % fitting, L0 should be the same as nlm0.LogLikelihood where 
            % nlm0 is created like this:
            %   nlm0 = fitnlm(ds,linfun0,rand);
            
            if isempty(model.Robust)
                % No Robust fitting. Get combined weights.
                w = get_CombinedWeights_r(model,false);
            else
                % Robust fitting. Get a vector of all ones.
                w = ones(length(model.ObservationInfo.Weights),1);
            end            

            % Subset of observations.
            subset = model.ObservationInfo.Subset;            
            % Subset of response.
            y = getResponse(model);
            y_r = y(subset);          
            % Subset of weights.
            w_r = w(subset);
            % Get weighted mean of y_r.
            mu0 = sum(y_r .* w_r)/sum(w_r);
            % Get weighted standard deviation around mu0.
            n_r = length(w_r);
            if ( n_r > 1 )
                mse0 = sum( w_r .* (y_r - mu0).^2 )/(n_r-1);
            else
                mse0 = sum( w_r .* (y_r - mu0).^2 )/(n_r);
            end
            L0 = -( (n_r-1) + n_r*log(2*pi) + sum(log(mse0./w_r)) )/2;
     
        end
        
        % --------------------------------------------------------------------
        function r = get_residuals(model,type)
            if nargin < 2 % return all residual types in a table array
                Raw = get_residuals(model,'raw');
                Pearson = get_residuals(model,'pearson');
                Studentized = get_residuals(model,'studentized');
                Standardized = get_residuals(model,'standardized');
                r = table(Raw,Pearson,Studentized,Standardized, ...
                            'RowNames',model.ObservationNames);
            else % return a single type of residual
                raw = getResponse(model) - predict(model);
                switch lower(type)
                case 'raw'
                    r = raw;
                case 'pearson'
                    r = raw ./ model.RMSE;
                case 'studentized' % "externally studentized", using Delete 1 Variances
                    h = model.Leverage;
                    s2_i = get_S2_i(model);
                    r = raw ./ sqrt(s2_i .* (1-h));
                case 'standardized' % "internally studentized", using MSE
                    r = raw ./ (sqrt(model.MSE * (1-model.Leverage)));
                otherwise
                    error(message('stats:NonLinearModel:BadResidualType', type));
                end
                r(~model.ObservationInfo.Subset) = NaN;
            end
        end
        
        function r = get_weighted_residuals(model,type)            
            % Get unweighted residuals.
            if ( nargin < 2 )
                r = get_residuals(model);
            else
                r = get_residuals(model,type);
            end
            
            % Get effective weights.
            if isempty(model.Robust)
                % No Robust fitting. Get combined weights.
                w = get_CombinedWeights_r(model,false);
            else
                % Robust fitting. Get a vector of all ones.
                w = ones(length(model.ObservationInfo.Weights),1);
            end

            % Multiply each column of r by sqrt(w).
            if isa(r,'dataset')
                r = dataset2table(r);
            end
            if isa(r,'table')
                rnames = r.Properties.VariableNames;
                r = varfun(@(xx) xx .* sqrt(w), r);
                r.Properties.VariableNames = rnames;                
            else
                r = sqrt(w) .* r;
            end   
        end
        
        
    end % protected methods
    
    methods(Static, Access='public', Hidden)
        function model = fit(X,varargin) % [X, y | DS], modelDef, initialCoefs, ...
% Not intended to be called directly. Use FITNLM to fit a NonLinearModel.
%
%   See also FITNLM.

            [X,y,haveDataset,otherArgs] = NonLinearModel.handleDataArgs(X,varargin{:});
            
            % Model definition and initial coefs must be given.
            if length(otherArgs) < 2
                error(message('stats:NonLinearModel:TooFewArguments'));
            end
            modelDef = otherArgs{1};
            initialCoefs = otherArgs{2};
            ncoefs = numel(initialCoefs);
            otherArgs(1:2) = [];
            
            % VarNames are optional names for the X matrix and y vector.  A
            % dataset/table defines its own list of names, so this is not accepted
            % with a dataset/table.
            
            paramNames = {'CoefficientNames' 'PredictorVars' 'ResponseVar' 'Weights' 'Exclude' ...
                          'VarNames' 'Options' 'ErrorModel' 'ErrorParameters'};
            paramDflts = {[] [] [] [] [] [] [] []};
            [coefNames,predictorVars,responseVar,weights, ...
             exclude,varNames,options,errormodel,errorparam,supplied] = ...
                internal.stats.parseArgs(paramNames, paramDflts, otherArgs{:});
            
            model = NonLinearModel();
            model.Iterative.InitialCoefs = initialCoefs;
            model.Iterative.IterOpts = options;
            if ~isempty(options)
                options = statset(statset('nlinfit'),options);
            end
            model.Robust = classreg.regr.FitObject.checkRobust(options);
            
            if ~haveDataset
                nx = size(X,2);
                [varNames,predictorVars,responseVar] = ...
                        classreg.regr.FitObject.getVarNames(varNames,predictorVars,responseVar,nx);
            end
            
            model.Formula = NonLinearModel.createFormula(supplied,modelDef,X,ncoefs,coefNames, ...
                predictorVars,responseVar,varNames,haveDataset);
            
           
            % If weights is a function handle, we will interpret it as an
            % error model. We will set model.ObservationInfo.Weights = all
            % ones but use GLS fitting in nlinfit.
            if ~isempty(weights) && isa(weights,'function_handle')
                % Make weights = [] in the call to assignData.
                model = assignData(model,X,y,[],[],model.Formula.VariableNames,exclude);                
                % Check weights. After data assignment, we can use getResponse(model).
                weights = NonLinearModel.checkWeights(weights, getResponse(model));
                % Store weights supplied as function handle in property WeightFunctionHandle.
                model.WeightFunctionHandle = weights;
            else
                % Either weights is a vector or we have no weights.
                model = assignData(model,X,y,weights,[],model.Formula.VariableNames,exclude);
            end
    
            if length(model.Formula.CoefficientNames) ~= ncoefs
                error(message('stats:NonLinearModel:BadInitial', ncoefs, length( model.Formula.CoefficientNames )));
            end
            
            % Validate errormodel and errorparam and store in ErrorModelInfo property.
            [errormodel,errorparam] = NonLinearModel.ValidateErrorModelAndErrorParam(errormodel,errorparam);
            model.ErrorModelInfo.ErrorModel = errormodel;
            model.ErrorModelInfo.ErrorParam = errorparam;
            
            % doFit is a FitObject method that does the following:
            % (1) model = selectVariables(model);            % Set PredLocs, RespLoc and update VariableInfo.InModel
            % (2) model = selectObservations(model,exclude); % Update ObservationInfo.Missing, .Excluded and .Subset
            % (3) model = fitter(model);                     % Do the actual fitting.
            % (4) model = postFit(model);                    % Do post fitting.           
            model = doFit(model);

            model = updateVarRange(model); % omit excluded points from range
        end
    end % static methods
    
    methods(Static = true, Hidden=true)
        function model = loadobj(obj)
            % Take care of base class changes first
            obj = loadobj@classreg.regr.ParametricRegression(obj);
            
            % (1) Set ErrorModelInfo   
            emptyErrorModel = false;
            if isempty(obj.ErrorModelInfo.ErrorModel)                                
                % Loading old object.
                emptyErrorModel = true; 
               
                % Initialize struct S.
                S = struct('ErrorModel',[],'ErrorParameters',[],'ErrorVariance',[],...
                    'MSE',[],'ScheffeSimPred',[],'WeightFunction',false,...
                    'FixedWeights',false,'RobustWeightFunction',false);

                % The 'constant' variance model, maybe using robust regression.
                mse = obj.MSE;
                errorparam = sqrt(mse);
                S.ErrorModel = 'constant';
                S.ErrorParameters = errorparam;
                S.ErrorVariance = @(x) mse * ones(size(x,1),1);  
                S.MSE = mse;
                S.ScheffeSimPred = obj.NumCoefficients + 1;                   
                if ~all(obj.ObservationInfo.Weights == 1)
                    % Fixed weights.
                    S.FixedWeights = true;
                end
                if ~isempty(obj.Robust)
                    % Robust weight function. 
                    S.RobustWeightFunction = true;
                end
                obj.ErrorModelInfo = S;
            end
            
            % (2) Set getEstimableContrasts - assume full rank.
            if isempty(obj.getEstimableContrasts)
                % Assume all contrasts are estimable.
                obj.getEstimableContrasts = @(Cnew) true(size(Cnew,1),1);
            end
            
            % (3) Set RankJW - assume full rank.
            if isempty(obj.RankJW)
                % Set equal to the number of coefficients.
                obj.RankJW = obj.NumCoefficients;
            end
            
            % (4) Set WeightFunctionHandle
            % Can be empty, no need to set.            
                        
            % (5) Set SST0 
            if isempty(obj.SST0) || isnan(obj.SST0)
                obj = getSumOfSquaresAndLogLikelihood(obj);
            end

            [f,p,emptyNullModel,hasIntercept] = fTest(obj);                      
            % (6) Set FTestFval
            if isempty(obj.FTestFval) || isnan(obj.FTestFval)
                obj.FTestFval = f;
            end
            
            % (7) Set FTestPval
            if isempty(obj.FTestPval) || isnan(obj.FTestPval)
                obj.FTestPval = p;
            end
            
            % (8) Set EmptyNullModel 
            if isempty(obj.EmptyNullModel) || isnan(obj.EmptyNullModel)
                obj.EmptyNullModel = emptyNullModel;
            end
            
            % (9) Set HasIntercept
            if isempty(obj.HasIntercept) || isnan(obj.HasIntercept)
                obj.HasIntercept = hasIntercept;
            end
               
            % (10) Re-set obj.ErrorModelInfo.ScheffeSimPred if required
            if emptyErrorModel && obj.HasIntercept
                obj.ErrorModelInfo.ScheffeSimPred = obj.NumCoefficients;
            end
            
            % Return
            model = obj;            
            
        end
    end

    
    methods(Static, Access='protected')
        function [X,y,haveDataset,otherArgs] = handleDataArgs(X,y,varargin) % [X, y | DS], ...
            if isa(X,'dataset')
                X = dataset2table(X);
            end
            haveDataset = isa(X,'table');
            if haveDataset
                % Have a dataset/table array, so no y.
                if nargin > 1
                    % Put the second arg in with the rest
                    otherArgs = [{y} varargin];
                else
                    otherArgs = {};
                end
                y = []; % just put anything in y
            elseif nargin < 2
                error(message('stats:NonLinearModel:MissingY'))
            else
                if isrow(X) && numel(X)==numel(y)
                    X = X(:);
                    y = y(:);
                end
                otherArgs = varargin;
            end
        end
        
        % ------------------- Externally-defined methods ---------------------
        function formula = createFormula(supplied,modelDef,X,ncoefs,coefNames,predictorVars,responseVar,varNames,haveDataset)
            
            givenFun = isa(modelDef,'function_handle');
            givenString = ~givenFun && internal.stats.isString(modelDef);
            
            if isa(modelDef,'classreg.regr.NonLinearFormula')
                formula = modelDef;
                
            elseif givenString || givenFun
                if givenString
                    if supplied.PredictorVars || supplied.ResponseVar
                        error(message('stats:NonLinearModel:NamesWithFormula'));
                    end
                elseif givenFun && classreg.regr.NonLinearFormula.isOpaqueFun(modelDef)
                    % If the model definition is an opaque function, we need both coef and
                    % predictor variable names to create a NonLinearFormula.  In the
                    % remaining cases (anonymous function, formula string, or formula
                    % object), NonLinearFormula can (try to) find names in the model
                    % definition itself if none were given, and will reconcile any
                    % specified names with the model definition.
                    %
                    % Since we have an initial coef vector, we can construct default coef
                    % names if none were given.
                    if ~supplied.CoefficientNames
                        coefNames = strcat({'b'},num2str((1:ncoefs)'));
                    end
                end
                
                if haveDataset
                    % VarNames is intended to name the columns of separate X and y.  It's
                    % not needed or accepted for a dataset/table array.
                    if supplied.VarNames
                        error(message('stats:NonLinearModel:NamesWithDataset'));
                    end
                    varNames = X.Properties.VariableNames;
                    
                else % data given as X,y
                    % ResponseVar is intended to say which variable in a dataset/table array is
                    % the response.  It's not needed or accepted for separate X and y.
                    nvars = size(X,2) + 1;
                    if supplied.VarNames
                        if length(varNames) ~= nvars
                            error(message('stats:NonLinearModel:BadVarNamesLength'));
                        end
                    else
                        if givenString
                            % If given a formula string and no names, NonLinearFormula will
                            % assume vars are in alphabetical order, but we have to force the
                            % response to be last.  Since we don't know the response name, we
                            % have to ask NonLinearFormula first.
                            formula = classreg.regr.NonLinearFormula(modelDef,coefNames,[],[],[],ncoefs);
                            varNames = formula.VariableNames;
                            responseName = formula.ResponseName;
                            respLoc = find(strcmp(responseName,varNames));
                            varNames = varNames([1:(respLoc-1) (respLoc+1):end respLoc]);
                            if length(varNames) ~= nvars
                                error(message('stats:NonLinearModel:CannotDetermineNames'));
                            end
                        else % givenFun
                            if classreg.regr.NonLinearFormula.isOpaqueFun(modelDef)
                                % Construct default variable names.
                                varNames = [strcat({'X'},num2str((1:size(X,2))'))' {'y'}];
                            end
                        end
                    end
                end
                
                if givenString
                    formula = classreg.regr.NonLinearFormula(modelDef,coefNames,[],[],varNames,ncoefs);
                else % givenFun
                    if classreg.regr.NonLinearFormula.isOpaqueFun(modelDef)
                        if ~supplied.PredictorVars
                            if haveDataset && supplied.ResponseVar
                                predictorVars = varNames(~strcmp(responseVar,varNames));
                            else
                                predictorVars = varNames(1:(end-1));
                            end
                       end
                    end
                    formula = classreg.regr.NonLinearFormula(modelDef,coefNames,predictorVars,responseVar,varNames,ncoefs);
                end
                
            else
                error(message('stats:NonLinearModel:BadModelDef'));
            end
        end
        
        %==== Subfunction ValidateErrorModelAndErrorParam ====
        function [errormodel, errorparam] = ValidateErrorModelAndErrorParam(errormodel,errorparam)
        % This function ensures that:
        % (1) the supplied errormodel is on the list of allowed error models and
        % (2) the supplied errorparam is a sensible value for the selected errormodel.
            
        
            % Make sure errormodel is one from the supported list.
            if isempty(errormodel)
                errormodel = NonLinearModel.NAME_CONSTANT;
            elseif ~isempty(errormodel) && ~any(strcmpi(errormodel,NonLinearModel.SupportedErrorModels))
                error(message('stats:NonLinearModel:InvalidErrorModel'));
            end
            
            % Make sure errorparam is sensible.
            if numel(errorparam)>2 || ~isnumeric(errorparam)
                error(message('stats:nlinfit:BadErrorParam'))
            end                        
            switch lower(errormodel)
                case NonLinearModel.NAME_COMBINED
                    if isempty(errorparam)
                        errorparam = [1 1];
                    elseif numel(errorparam)~= 2
                        % For combined error model, errorparam should be a vector [a b].
                        error(message('stats:nlinfit:BadCombinedParam', errormodel));
                    end
                case NonLinearModel.NAME_PROPORTIONAL
                    % Only a should be specified.
                    if isempty(errorparam)
                        errorparam = 1;
                    elseif numel(errorparam)~=1
                        error(message('stats:nlinfit:BadErrorParam1', errormodel))
                    end
                case NonLinearModel.NAME_CONSTANT
                    % Only b should be specified.
                    if isempty(errorparam)
                        errorparam = 1;
                    elseif numel(errorparam)~=1
                        error(message('stats:nlinfit:BadErrorParam1', errormodel))
                    end
            end
            
        end % End of ValidateErrorModelAndErrorParam.
            
        %==== Subfunction checkWeights ====
        function weights = checkWeights(weights, yfit)
        % If weights is a function handle, make sure it works as required. Whether
        % weights vector is given directly or through a function handle, make sure 
        % it is numeric, real, positive and of the same length as the number of 
        % predictions.

            if isa(weights, 'function_handle') || (isnumeric(weights) && isvector(weights))
                % If weights are set.

                if isa(weights, 'function_handle')
                    % function handle.
                    try
                        wVec = weights(yfit);
                    catch ME
                        if isa(weights, 'inline')
                            m = message('stats:nlinfit:InlineWeightFunctionError');
                            throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
                        elseif strcmp('MATLAB:UndefinedFunction', ME.identifier) ...
                                && ~isempty(strfind(ME.message, func2str(weights)))
                            error(message('stats:nlinfit:WeightFunctionNotFound',func2str(weights)));
                        else
                            m = message('stats:nlinfit:WeightFunctionError',func2str(weights));
                            throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
                        end
                    end
                else
                    % fixed weights.
                    wVec = weights;
                end

                % wVec should be real positive vector of the same size as yfit.
                if ~isequal(size(wVec), size(yfit)) || ~isnumeric(wVec) || ~isreal(wVec) || ~isvector(wVec) || any( wVec(:) <= 0 ) 
                    error(message('stats:nlinfit:InvalidWeights'));
                end

            else
               % Invalid weights.
               error(message('stats:nlinfit:InvalidWeights'));
            end

        end % End of checkWeights.
        
        %==== Subfunction validateScheffeParamUsingJacobianPred ====
        function sch = validateScheffeParamUsingJacobianPred(sch, rankJ, delta, estimable, errorVar)
        % Check if Jacobians used during model fitting and prediction have the same structure.

            if sch ~= (rankJ+1)                
                % Get only estimable rows in delta.
                delta_est = delta(estimable,:);
                % Since EVFit = ones(size(EVPred)), Unweighted Jacobian = Weighted Jacobian.
                EVPred = errorVar(estimable); EVFit = ones(size(EVPred));                
                [schEstDelta,rankEstDelta] = internal.stats.getscheffeparam('UnWeightedJacobian',delta_est,'Intopt','observation','ErrorVarianceFit',EVFit,'ErrorVariancePred',EVPred);
                if ( schEstDelta ~= rankEstDelta )
                   % Prediction Jacobian has a different structure. If Jpred = Vr*Sigmar*Ur' 
                   % then Vr*Vr'*dpred \neq dpred where dpred(i) = sqrt(Gpred_ii). Reset sch
                   % to (rankJ+1).
                   sch = (rankJ+1); 
                end        
            end

        end % End of validateScheffeParamUsingJacobianPred.        
        
        %==== Subfunction applyLogTransformation ====
        function [y, model] = applyLogTransformation(y, model)
            % Exponential error model. Linearize the model as
            %   y = f*exp(a*e), or log(y) = log(f) + a*e

            if ~all(y>0)
                error(message('stats:nlinfit:PositiveYRequired'));
            else
                y = log(max(y,realmin));
            end

            model = @(phi,X) log(max(model(phi,X),realmin));
        end % End of applyLogTransformation.
        
    end % static protected methods        
    
end
