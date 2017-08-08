classdef (Sealed = true) GeneralizedLinearMixedModel < classreg.regr.LinearLikeMixedModel
%GeneralizedLinearMixedModel Fitted generalized linear mixed effects model.
%   GLME = FITGLME(...) fits a generalized linear mixed effects model to
%   data. The fitted model GLME is a GeneralizedLinearMixedModel object
%   modeling the response variable from one of several supported
%   distributions. A transformation of the conditional mean of the response
%   is modelled as a linear function of fixed and random effect predictors.
%   Several choices for the transformation are available via the supported
%   link functions.
%
%   GeneralizedLinearMixedModel methods:
%       coefCI - Coefficient confidence intervals
%       coefTest - Linear hypothesis test on coefficients
%       predict - Compute predicted values given predictor values
%       random - Generate random response values given predictor values
%       plotResiduals - Plot of residuals
%       designMatrix - Fixed and random effects design matrices
%       fixedEffects - Stats on fixed effects
%       randomEffects - Stats on random effects
%       covarianceParameters - Extract covariance parameters
%       fitted - Fitted response
%       residuals - Various types of residuals
%       response - Response used to fit the model
%       compare - Compare fitted models
%       anova - Marginal tests for fixed effect terms
%       disp - Display a fitted model
%       refit - Fit the model to new response vector
%
%   GeneralizedLinearMixedModel properties:
%       FitMethod - Method used for fitting
%       Distribution - Distribution of the response
%       Link - Link relating the distribution parameters to the predictors
%       Dispersion - Theoretical or estimated dispersion parameter
%       DispersionEstimated - Flag indicating if Dispersion was estimated
%       Formula - Representation of the model used in this fit
%       LogLikelihood - Log of likelihood function at coefficient estimates
%       DFE - Degrees of freedom for residuals
%       SSE - Error sum of squares
%       SST - Total sum of squares
%       SSR - Regression sum of squares
%       CoefficientCovariance - Covariance matrix for coefficient estimates
%       CoefficientNames - Coefficient names
%       NumCoefficients - Number of coefficients
%       NumEstimatedCoefficients - Number of estimated coefficients
%       Coefficients - Coefficients and related statistics
%       Rsquared - R-squared and adjusted R-squared
%       ModelCriterion - AIC and other model criteria
%       VariableInfo - Information about variables used in the fit
%       ObservationInfo - Information about observations used in the fit
%       Variables - Table of variables used in fit
%       NumVariables - Number of variables used in fit
%       VariableNames - Names of variables used in fit
%       NumPredictors - Number of predictors
%       PredictorNames - Names of predictors
%       ResponseName - Name of response
%       NumObservations - Number of observations in the fit
%       ObservationNames - Names of observations in the fit
%
%   See also FITGLME, FITLME, LinearMixedModel.    
     
%   Copyright 2013-2014 The MathWorks, Inc.

% Properties from Predictor.
properties(Dependent,GetAccess='public',SetAccess='protected',Hidden=true)    
%Fitted - Fitted (predicted) values.
%   The Fitted property is a vector of estimated conditional means that
%   include contributions from both fixed and random effects and are on the
%   scale of the response variable.
%
%   The fitted values are computed using the predictor values used to fit
%   the model. Use the predict method to compute predictions for other
%   predictor values and to compute confidence bounds for the predicted
%   values.
%
%   See also FITTED.
    Fitted
    
%Residuals - Residual values.
%   The Residuals property is a table array with the following columns:
%
%        'Raw'          The Raw conditional residuals
%        'Pearson'      The Pearson conditional residuals
%
%   The conditional residuals are based on fitted values that include
%   contributions from both fixed and random effects and are on the scale
%   of the response variable.
%
%   For more details on various residuals, see the residuals method.
%
%   See also RESIDUALS.
    Residuals    
end
  
% New public properties defined by GeneralizedLinearMixedModel.
properties(GetAccess=public, SetAccess='protected')
%FitMethod - Method used to fit the generalized linear mixed effects model.
%   The FitMethod property is one of the following: 
%
%                  'MPL' - Maximum pseudo likelihood
%                'REMPL' - Restricted maximum pseudo likelihood
%   'ApproximateLaplace' - Maximum likelihood using approximate Laplace
%                          approximation with fixed effects profiled out
%              'Laplace' - Maximum likelihood using Laplace approximation
%
%   See also FITGLME.
    FitMethod
    
    
%Distribution Response distribution name.
%   The Distribution property is a string containing the name of the
%   distribution of the response.
%
%   See also GeneralizedLinearMixedModel, Link, Dispersion, DispersionEstimated.
    Distribution    
    
%Link Link between the distribution parameters and predictor values.
%   The Link property is a structure providing the name and other
%   characteristics of the link function. The link is a function G that
%   links the distribution parameter MU to the linear predictor ETA as
%   follows: G(MU) = ETA. The structure has four fields:
%
%        Name                 Name of the link function.
%        Link                 The function that defines G.
%        Derivative           Derivative of G.
%        SecondDerivative     Second derivative of G.
%        Inverse              Inverse of G.
%
%   See also GeneralizedLinearMixedModel, Distribution.
    Link    
    
%Dispersion Parameter defining the conditional variance of the response.
%   For observation i, the conditional (given the random effects) variance
%   of the response y_i given the conditional mean mu_i and dispersion
%   parameter phi in a generalized linear mixed effects (GLME) model is
%   given by:
%
%           Var(y_i | mu_i, phi) = (phi/w_i) * v(mu_i),      phi > 0
%
%   where w_i is the i th observation weight and v is the variance function
%   for the specified conditional distribution of the response. The
%   Dispersion property contains an estimate of phi for the specified GLME
%   model. The value of Dispersion depends on the specified conditional
%   distribution of the response. For binomial and Poisson distributions,
%   the theoretical value of Dispersion is equal to phi = 1.0.
%
%   If 'FitMethod' is 'MPL' or 'REMPL' and the 'DispersionFlag' name/value
%   pair in fitglme is true then a dispersion parameter is estimated from
%   data even for binomial and Poisson distributions. If 'FitMethod' is
%   equal to 'ApproximateLaplace' or 'Laplace' then 'DispersionFlag'
%   name/value pair in fitglme does not apply and the dispersion parameter
%   is always fixed at 1.0 for binomial and Poisson distributions. For all
%   other distributions, Dispersion is always estimated from data.
%
%   See also GeneralizedLinearMixedModel, DispersionEstimated.
    Dispersion        
        
%DispersionEstimated A flag indicating if dispersion was estimated.
%   For the binomial and Poisson distributions the dispersion parameter has
%   a known theoretical value of 1.0. If FitMethod is 'ApproximateLaplace'
%   or 'Laplace', the dispersion parameter is fixed at its theoretical
%   value of 1.0 for binomial and Poisson distributions.
%
%   If FitMethod is 'MPL' or 'REMPL', the dispersion parameter is estimated
%   even for binomial and Poisson distributions if the 'DispersionFlag'
%   name/value pair in fitglme is true and fixed at its theoretical value
%   of 1.0 if the 'DispersionFlag' name/value pair in fitglme is false. 
%
%   For distributions other than binomial and Poisson and for any
%   'FitMethod', the dispersion parameter is always estimated.
%
%   See also GeneralizedLinearMixedModel, Dispersion.
    DispersionEstimated                
end

% New public hidden properties defined by GeneralizedLinearMixedModel.
properties(GetAccess='public', SetAccess='protected',Hidden=true)    
%Response - Vector of response values used to fit the model
    Response    
end

% Internal constants.
properties(Constant=true,Hidden=true)        
    
    AllowedFitMethods = {'mpl','rempl','approximatelaplace','laplace'};        
        
    AllowedOptimizers = {'fminunc','quasinewton','fminsearch'};

    AllowedStartMethods = {'random','default'};         

    AllowedDFMethods = {'none','residual'};
        
    AllowedResidualTypes = {'raw','pearson'};
    
    AllowedDistributions = {'normal','gaussian','binomial','poisson','gamma','inversegaussian','inverse gaussian'};    

    AllowedEBMethods = {'auto','default','linesearchnewton','linesearchmodifiednewton','trustregion2d','fsolve'}; 
    
    AllowedCovarianceMethods = {'conditional','jointhessian'};
                            
end

% Private properties.
properties(Access={?classreg.regr.LinearLikeMixedModel})
%   The following variables only include observations used in the fit 
%   corresponding to ObservationInfo.Subset.
    
%BinomialSize - Vector of integers containing the number of trials using 
%   which the proportions in y were computed. Applies only when modeling a
%   binomial distribution.
        BinomialSize    
        
%Offset - Offset vector for the linear predictor of the GLME model.
        Offset
        
%DispersionFlag - Supplied value of 'DispersionFlag' property.    
        DispersionFlag 
    
%CovariancePattern - A cell array of covariance patterns.
%   If the GLME model has R grouping variables, CovariancePattern is a cell
%   array of length R containing the covariance pattern used for each 
%   grouping variable. 
        CovariancePattern
    
%DummyVarCoding - A string containing the dummy variable coding used for 
%   creating design matrices.
        DummyVarCoding
    
%Optimizer - A string containing the optimizer to use for optimization.
        Optimizer
    
%OptimizerOptions - A structure or object specifying the optimization 
%   options.
        OptimizerOptions
    
%StartMethod - A string specifying how to initialize parameters for 
%   optimization.
        StartMethod
    
%CheckHessian - A logical scalar indicating whether Hessian checks should 
%   be performed after fitting the model.    
        CheckHessian

%PLIterations - A positive integer specifying the maximum number of pseudo 
%   likelihood (PL) iterations. Applies only if 'FitMethod' is 'MPL' or 
%   'REMPL'.
        PLIterations    
        
%PLTolerance - A numeric, real scalar indicating the relative convergence
%   tolerance for terminating PL iterations.
        PLTolerance   
        
%MuStart - A vector providing a starting value for the conditional mean of 
%   y given b to initialize PL iterations. Legal values of elements in 
%   MuStart are as follows:
%
%                 Distribution          Legal values
%                 'normal'              (-Inf,Inf)
%                 'binomial'            (0,1)
%                 'poisson'             (0,Inf)
%                 'gamma'               (0,Inf)
%                 'inverse gaussian'    (0,Inf)        
%
%   If MuStart is empty, a default initialization for MuStart is used.
        MuStart
    
%InitPLIterations - Initial number of PL iterations used to initialize ML 
%   based methods like 'Laplace' and 'ApproximateLaplace'.
        InitPLIterations
        
%EBMethod - A string containing the optimization method to use for
%   estimating the empirical Bayes (EB) predictors of random effects. 
%   Choices are 'Default', 'LineSearchNewton', 'TrustRegion2D' and 
%   'fsolve'. The 'Default' 'EBMethod' is very similar to 
%   'LineSearchNewton' but uses a different convergence criterion and does
%   not display iterative progress. 'Default' and 'LineSearchNewton' may
%   fail for non-canonical link functions. In these cases, 'TrustRegion2D'
%   or 'fsolve' is the recommended 'EBMethod'. 'fsolve' can be used only if
%   we detect Optimization Toolbox.
        EBMethod     
        
%EBOptions - A structure created using statset to control the EB
%   optimization. EBOptions has the following fields:
%
%                 'TolFun'      Relative tolerance on the gradient norm.
%                               Default is 1e-6.
%                 'TolX'        Absolute tolerance on the step size.
%                               Default is 1e-8.
%                 'MaxIter'     Maximum number of iterations. Default is
%                               100.
%                 'Display'     'off', 'iter' or 'final'. Default is 'off'.
%
%   For 'EBMethod' equal to 'Default', 'TolFun' is the relative tolerance
%   on the linear predictor of the model and the option 'Display' does not
%   apply. If 'EBMethod' is 'fsolve', this should be an object created by
%   optimoptions('fsolve').
        EBOptions                         

%CovarianceMethod - A string indicating the method used to compute the 
%   covariance matrix of estimated parameters. Choices are 'Conditional'
%   and 'JointHessian'. The 'Conditional' method is fast and computes an
%   approximate covariance of fixed effects conditional on the estimated
%   covariance parameters equal to the true values. The 'JointHessian'
%   method computes the joint Hessian of the Laplacian log likelihood with
%   respect to all model parameters and uses this to compute the covariance
%   of fixed effects and covariance parameters. Default is 'Conditional'.
%   For the 'Conditional' method, the covariance of covariance parameters
%   is not computed.
        CovarianceMethod        
        
%UseSequentialFitting - Either true or false. If false, all ML methods are 
%                       initialized using 1 or more PL iterations. If true,
%                       the initial values from PL iterations are refined
%                       with 'ApproximateLaplace' for 'Laplace' based
%                       fitting.
        UseSequentialFitting                       
    
%CovarianceTable - A cell array of tables returned by the method
%   covarianceParameters. We store this info here so that we don't have to
%   recompute this in disp.
        CovarianceTable
        
%ShowPLOptimizerDisplay - Logical flag indicating whether we want to show 
%   the iterative display from PL optimizer.    
        ShowPLOptimizerDisplay = false;        
end

% Abstract public methods from FitObject.
methods(Access='public', Hidden=true)
    
    function t = title(model)
        strLHS = model.ResponseName;
        strFunArgs = internal.stats.strCollapse(model.Formula.PredictorNames,',');
        t = sprintf( '%s = glme(%s)',strLHS,strFunArgs);
    end % end of title.        
    
    function val = feval(model,varargin) %#ok<INUSD>
        warning(message('stats:GeneralizedLinearMixedModel:NoFevalMethod'));
        val = [];
    end % end of feval.
        
end

% Abstract public methods from FitObject.
methods(Access='public')
        
    function disp(model)
%DISP Display a GeneralizedLinearMixedModel.
%   DISP(GLME) displays the GeneralizedLinearMixedModel object GLME.
%
%   See also GeneralizedLinearMixedModel, FITGLME.

        % 1. Can't display an empty model.
        if isempty(model.ObservationInfo)
            displayFormula(model);
            error(message('stats:GeneralizedLinearMixedModel:NoConstructor'));
        end
        
        % 2. Display headline.
        displayHeadLine(model);
        
        % 3. Model information.
        displayModelInfo(model);
        
        % 4. Display the formula.
        displayFormula(model);

        % 5. Model fit stats.
        displayModelFitStats(model);        
            
        % 6. Fixed effect stats.
        displayFixedStats(model)            
            
        % 7. Random effects covariance parameter stats.
        displayCovarianceStats(model);               
        
    end % end of disp.
        
end

% Abstract public methods from Predictor.
methods(Access='public')

    function [mupred,muci,df] = predict(model,varargin)
%PREDICT Compute predicted values given predictor values.
%   MUPRED = PREDICT(GLME) computes a vector MUPRED of predictions of the
%   conditional mean of the response given the random effects at the
%   original predictors used to create GLME. MUPRED includes contributions
%   from both estimated fixed effects and approximate empirical Bayes
%   predictors (EBPs) of the random effects.
%
%   MUPRED = PREDICT(GLME,T) uses predictor variables from the table T to
%   predict the conditional mean of the response. T must contain all of the
%   predictor variables used to create GLME.
%
%   [MUPRED,MUCI] = PREDICT(...) also returns the two-column matrix MUCI
%   containing 95% pointwise confidence intervals for the predicted values.
%   The lower limits of the bounds are in column 1, and the upper limits
%   are in column 2.
%
%   [MUPRED,MUCI,DF] = PREDICT(...) also returns the degrees of freedom DF.
%   If the 'Simultaneous' name-value pair is 'false', then DF is a vector
%   containing the degrees of freedom values used in computing the
%   confidence intervals. If 'Simultaneous' is 'true', then DF is a scalar
%   containing the degrees of freedom value used for computing the
%   simultaneous confidence intervals.
%  
%   [...] = PREDICT(...,PARAM1,VALUE1,...) specifies one or more of the
%   following name/value pairs:
%
%     'Conditional'     Either true or false. If false then returns
%                       marginal predictions which include contribution
%                       only from the estimated fixed effects. When making
%                       predictions, if a particular grouping variable has
%                       new levels (ones that were not in the original
%                       data), then the random effects for the grouping
%                       variable do not make any contribution to the
%                       'Conditional' prediction at observations where the
%                       grouping variable has new levels. Default is true.
%
%    'Simultaneous'     Either true for simultaneous bounds, or false for 
%                       non-simultaneous bounds. Default is false.
%
%       'DFMethod'      Specifies the method to use for computing the
%                       approximate denominator degrees of freedom (DF)
%                       when computing the confidence intervals. Options
%                       are 'Residual' and 'None'. If 'DFMethod' is
%                       'Residual', the DF values are assumed to be
%                       constant and equal to (N-P) where N is the number
%                       of observations and P is the number of fixed
%                       effects. If 'DFMethod' is 'None', then all DF
%                       values are set to infinity. Default is 'Residual'.
%  
%          'Alpha'      A value between 0 and 1 to specify the confidence
%                       level as 100(1-ALPHA)%.  Default is 0.05 for 95%
%                       confidence.
%
%           'Offset'    Specifies the offset vector of the model. Must be
%                       a vector of length M where M is the number of rows 
%                       in T. Default is zeros(M,1).
%
%   predict computes the confidence intervals using the conditional mean
%   squared error of prediction (CMSEP) approach of Booth and Hobert (1998)
%   conditional on the estimated covariance parameters and the observed
%   response. An alternative interpretation of the confidence intervals
%   from predict is that they are approximate Bayesian credible intervals
%   conditional on the estimated covariance parameters and the observed
%   response.
%
%   The marginal predictions are not the predictions of the marginal mean 
%   of the response. Marginal predictions are computed by substituting a
%   vector of zeros in place of the EBP vector of the random effects in the
%   prediction for the conditional mean of the response given the random
%   effects.
%
%   Example: Model gas mileage as a function of car weight, with a random
%            effect due to model year.
%      load carsmall
%      T = table(MPG,Weight,Model_Year);
%      glme = fitglme(T,'MPG ~ Weight + (1|Model_Year)')
%
%      % Plot predicted values conditional on each year.
%      gscatter(T.Weight,T.MPG,T.Model_Year)
%      T2 = table();
%      T2.Weight = linspace(1500,5000)';
%      T2.Model_Year = repmat(70,100,1);
%      line(T2.Weight,predict(glme,T2),'color','r');
%      T2.Model_Year(:) = 76;
%      line(T2.Weight,predict(glme,T2),'color','g');
%      T2.Model_Year(:) = 82;
%      line(T2.Weight,predict(glme,T2),'color','b');
%
%   See also FITTED, RANDOM.

%   References:  
%   (1) James G. Booth and James P. Hobert (1998), Standard Errors of
%   Prediction in Generalized Linear Mixed Models, Journal of the American
%   Statistical Association, 93(441), 262-272.

        if nargin < 2 || internal.stats.isString(varargin{1})            
            % 1. If no new design points are given, get conditional
            % predictions from the fitted model at the original design
            % points. Called like: 
            % mupred = predict(glme) or 
            % mupred = predict(glme,'Conditional',true,...).
            haveDataset = true;
            ds = model.Variables;            
            otherArgs = varargin;     
        else                        
            % 2. Get supplied dataset and optional args. Called like:
            % mupred = predict(glme,T,'Conditional',true,...).
            [haveDataset,ds,~,~,~,otherArgs] ...
                = model.handleDatasetOrMatrixInput(varargin{:});
        end
        assertThat(haveDataset,'stats:LinearMixedModel:MustBeDataset','T')
        M = size(ds,1);
        
        % 3. Parse and validate optional input args.
            % 3.1 Default parameter values.
            dfltConditional  = true;
            dfltSimultaneous = false;
            dfltDFMethod     = 'Residual';
            dfltAlpha        = 0.05;
            dfltOffset       = zeros(M,1);
            % 3.2 Optional parameter names and their default values.
            paramNames =   {'Conditional',   'Simultaneous',   'DFMethod',   'Alpha',   'Offset'};
            paramDflts = {dfltConditional, dfltSimultaneous, dfltDFMethod, dfltAlpha, dfltOffset};
            % 3.3 Parse optional parameter name/value pairs.
            [conditional,simultaneous,dfmethod,alpha,offset] ...
                = internal.stats.parseArgs(paramNames,paramDflts,otherArgs{:});   
            % 3.4 Validate optional parameter values.
            conditional  = model.validateConditional(conditional);
            simultaneous = model.validateSimultaneous(simultaneous);
            dfmethod     = model.validateDFMethod(dfmethod);
            alpha        = model.validateAlpha(alpha);
            offset       = model.validateOffset(offset,M);
  
        % 4. Extract size info from model.
          q = model.RandomInfo.q;
          R = model.GroupingInfo.R;
        lev = model.GroupingInfo.lev;
            
        % 5. Validate ds.        
        predNames = model.PredictorNames;
        varNames = model.Variables.Properties.VariableNames;
        [~,predLocs] = ismember(predNames,varNames);
        dsref = model.Variables(:,predLocs);
        ds = model.validateDataset(ds,'DS',dsref);        
        
        % 6. Get X and Z.
            % 6.1 Get X.
            finfo = extractFixedInfo(model,ds);    
                X = finfo.X;
            % 6.2 Get Z.
            rinfo = extractRandomInfo(model,ds);
                Z = rinfo.Z; 
                
        % 7. Get cell arrays Gid and GidLevelNames for the supplied
        % observations.
        ginfo = extractGroupingInfo(model,ds);
        Gid = ginfo.Gid;
        GidLevelNames = ginfo.GidLevelNames;       
        
        % 8. Modify Gid to newGid such that in newGid{i}, ID j is mapped
        % to model.GroupingInfo.GidLevelNames{i}{j}. In other words, if we
        % associated ID 5 with level name 'Green', we must associate ID 5
        % with level 'Green' in the supplied observations as well. New
        % levels in supplied data, not previously seen will be assigned a
        % group ID of NaN.
        newGid = cell(R,1);
        for k = 1:R           
            newGid{k} = model.reorderGroupIDs(Gid{k},...
                GidLevelNames{k},model.GroupingInfo.GidLevelNames{k});
        end
        
        % 9. X already contains the fixed effects design matrix. Get the
        % full sparse random effects design matrix Zs. This matrix will
        % have sum_{i=1 to R} (q(i)*lev(i)) columns and M rows.
        Zs = model.makeSparseZ(Z,q,lev,newGid,M);
  
        % 10. Get ypred and yci.        
            % 10.1 Prediction options.            
            wantConditional = conditional;
            if simultaneous == true
                wantPointwise = false;
            else
                wantPointwise = true;
            end                
            % 10.2 Call predict on model.slme.
            hasIntercept = model.Formula.FELinearFormula.HasIntercept;
            args = {X,Zs,alpha,dfmethod,...
                wantConditional,wantPointwise,offset,hasIntercept};
            switch nargout
                case {0,1}
                    mupred           = predict(model.slme,args{:});
                case 2
                    [mupred,muci]    = predict(model.slme,args{:});
                case 3
                    [mupred,muci,df] = predict(model.slme,args{:});
            end
            %[mupred,muci,df] = predict(model.slme,X,Zs,alpha,dfmethod,...
            %    wantConditional,wantPointwise,offset,hasIntercept);
                                    
    end % end of predict.
    
    function ynew = random(model,varargin)
%RANDOM Generate random response values given predictor values.
%   YNEW = RANDOM(GLME) generates a vector YNEW of random values from the
%   fitted generalized linear mixed effects model GLME at the original
%   design points. YNEW is created by first generating the random effects
%   vector from its fitted prior distribution and then generating YNEW from
%   its fitted conditional distribution given the random effects. The
%   effect of supplied observation weights when creating the fit (if any)
%   is taken into account.
%
%   YNEW = RANDOM(GLME,TNEW) uses predictor variables from the table TNEW,
%   which must contain all of the predictor variables used to create GLME.
%
%   YNEW = RANDOM(...,PARAM1,VALUE1,...) accepts optional name/value pairs:
%
%           'Name'       'Value'
%          'Weights'     Vector of M non-negative weights, where M is the
%                        number of rows in TNEW. Default is ones(M,1). For
%                        binomial and Poisson distributions, 'Weights' must
%                        be a vector of positive integers.
%  
%     'BinomialSize'     Vector of length M, where M is the number of rows
%                        in TNEW. Default is ones(M,1). Applies only to the
%                        binomial distribution and specifies the number of
%                        binomial trials when generating the random
%                        response values. Elements of 'BinomialSize' must
%                        be positive integers.
%
%           'Offset'     Specifies the offset vector of the model. Must be
%                        a vector of length M where M is the number of rows 
%                        in TNEW. Default is zeros(M,1).
%
%   Example: Fit a model. Simulate new random values for the first
%            observation.
%      load carsmall
%      T = table(MPG,Weight,Model_Year);
%      glme = fitglme(T,'MPG ~ Weight + (1|Model_Year)');
%      newT = T(ones(10000,1),:);
%      r = random(glme,newT);
%      hist(r)
%
%      % These random values share a common new random effect due to
%      % Model_Year, but their variance is comparable to the model value.
%      v1 = var(r)
%      v2 = glme.Dispersion
%
%   See also PREDICT, FITTED.

        % 1. If no new design points are given, generate data from the
        % fitted model at the original design points. Include the effect of
        % originally supplied observation weights.
        if nargin < 2                        
            % Called like: ynew = random(model).
            wp      = model.slme.PriorWeights;
            delta   = model.slme.Offset;
            ntrials = model.slme.BinomialSize;            
            ysim    = random(model.slme,[],model.slme.X,model.slme.Z,delta,wp,ntrials);                        
            
            subset       = model.ObservationInfo.Subset;
            ynew         = NaN(length(subset),1);       
            ynew(subset) = ysim;                        
            return;            
        end

        % 2. Should we use the original dataset or a new one?
        if internal.stats.isString(varargin{1})
            % 2.1 Called like: ynew = random(model,'Weights',w,'Offset',offset,'BinomialSize',binomsize).
            haveDataset = true;
            ds = model.Variables;            
            otherArgs = varargin;
        else
            % 2.2 Called like: ynew = random(model,Tnew,'Weights',w,'Offset',offset,'BinomialSize',binomsize).
            [haveDataset,ds,~,~,~,otherArgs] = model.handleDatasetOrMatrixInput(varargin{:});
        end                         
        assertThat(haveDataset,'stats:LinearMixedModel:MustBeDataset','TNEW');  
        M = size(ds,1);
        
        % 3. Validate optional name/value pairs.        
            % 3.1 Default parameter values.
            dfltWeights      = ones(M,1);
            dfltBinomialSize = ones(M,1);
            dfltOffset       = zeros(M,1);
            % 3.2 Optional parameter names and their default values.
            paramNames =   {'Weights',   'BinomialSize',   'Offset'};
            paramDflts = {dfltWeights, dfltBinomialSize, dfltOffset};
            % 3.3 Parse optional parameter name/value pairs.
            [weights,binomsize,offset] = internal.stats.parseArgs(paramNames,paramDflts,otherArgs{:});
            % 3.4 Validate optional parameter values.
            weights   = model.validateWeights(weights,M,model.Distribution);
            binomsize = model.validateBinomialSize(binomsize,M);
            offset    = model.validateOffset(offset,M);          
        
        % 4. Validate ds.        
        predNames = model.PredictorNames;
        varNames = model.Variables.Properties.VariableNames;
        [~,predLocs] = ismember(predNames,varNames);
        dsref = model.Variables(:,predLocs);
        ds = model.validateDataset(ds,'DS',dsref);
                                            
        % 5. Extract size info from model.
        q = model.RandomInfo.q;                                                                               
        
        % 6. Get X and Z from ds.
            % 6.1 Get X.
            finfo = extractFixedInfo(model,ds);    
                X = finfo.X;
            % 6.2 Get Z.
            rinfo = extractRandomInfo(model,ds);
                Z = rinfo.Z;  
                
        % 7. Get a R-by-1 cell vector Gid such that Gid{i} is the group
        % IDs for grouping variable i. Also get a R-by-1 vector lev such
        % that lev(i) is the number of levels of grouping variable i in the
        % input data. Extract grouping variable info from ds. In ginfo,
        % GidLevelNames may contain levels not previously seen during model
        % fitting. That's okay, every level gets its own random effects
        % vector.
        ginfo = extractGroupingInfo(model,ds);
        Gid = ginfo.Gid;
        lev = ginfo.lev;        
            
        % 8. X already contains the fixed effects design matrix. Get the
        % full sparse random effects design matrix Zs. This matrix will
        % have sum_{i=1 to R} (q(i)*lev(i)) columns and M rows.
        Zs = model.makeSparseZ(Z,q,lev,Gid,M);        
        
        % 9. Call random on the standard object. We want the random effects
        % vector to be of size sum_{i=1 to R} (q(i)*lev(i)) - so pass in
        % lev as the last input argument.
        ynew = random(model.slme,[],X,Zs,offset,weights,binomsize,lev);
        
    end % end of random.
    
end

% Abstract public methods from ParametricRegression.
methods(Access='public', Hidden=true)
    
    function v = varianceParam(model) % a consistent name for the dispersion/MSE/whatever
        v = model.Dispersion;
    end
    
end

% Public methods from ParametricRegression.
methods(Access='public')
   
    function [feci,reci] = coefCI(model,varargin)
%coefCI Confidence intervals for coefficients.
%   FECI = coefCI(GLME) computes 95% confidence intervals for the fixed
%   effects parameters in the generalized linear mixed effects model GLME.
%   The output FECI is a P-by-2 matrix where P is the number of fixed
%   effects parameters in GLME. Rows of FECI from top to bottom correspond
%   respectively to the P-by-1 fixed effects vector BETA displayed from top
%   to bottom in the tabular display from the fixedEffects method. Column 1
%   of FECI displays lower confidence limits and column 2 of FECI displays
%   upper confidence limits.
%
%   [FECI,RECI] = coefCI(GLME) also returns 95% confidence intervals for
%   the random effects parameters in GLME. The output RECI is a Q-by-2
%   matrix where Q is the total number of random effects parameters in
%   GLME. Rows of RECI from top to bottom correspond respectively to the
%   Q-by-1 random effects vector B displayed from top to bottom in the
%   tabular display from the randomEffects method. Column 1 of RECI
%   displays lower confidence limits and column 2 of RECI displays upper
%   confidence limits.
%
%   [FECI,RECI] = coefCI(GLME,'PARAM1','VALUE1',...) accepts optional
%   name/value pairs:
%      
%           'Name'     'Value'
%          'Alpha'     ALPHA, a number between 0 and 1. Computes 
%                      100*(1-ALPHA)% confidence intervals. Default is
%                      ALPHA=0.05 for 95% confidence intervals.
%  
%       'DFMethod'     Specifies the method to use for computing the
%                      approximate degrees of freedom DF for confidence
%                      interval calculation. Options are 'Residual' and
%                      'None'. If 'DFMethod' is 'Residual', the DF values
%                      are assumed to be constant and equal to (N-P) where
%                      N is the number of observations and P is the number
%                      of fixed effects. If 'DFMethod' is 'None', then all
%                      DF values are set to infinity. Default is
%                      'Residual'.
%
%   coefCI uses an approximation to the conditional mean squared error of
%   prediction (CMSEP) of the random effects to compute confidence
%   intervals. This accounts for the uncertainty in the estimation of
%   fixed effects but not the covariance parameters. An alternative
%   interpretation of RECI is that they are Bayesian credible intervals
%   using an approximation to the posterior distribution of random effects
%   given the estimated covariance parameters and the response.
%
%   Example: Fit a model with Weight as a fixed effect, and a random effect
%            due to Model_Year. Compute confidence intervals for the
%            intercept and the coefficient of Weight.
%      load carsmall
%      T = table(MPG,Weight,Model_Year);
%      glme = fitglme(T,'MPG ~ Weight + (1|Model_Year)','Distribution','Normal')
%      coefCI(glme)
%
%   See also coefTest, fixedEffects, randomEffects, covarianceParameters.
        
        switch nargout
            case {0,1}
                feci        = coefCI@classreg.regr.LinearLikeMixedModel(model,varargin{:});
            case 2
                [feci,reci] = coefCI@classreg.regr.LinearLikeMixedModel(model,varargin{:});
        end

    end % end of coefCI.

    function [P,F,DF1,DF2] = coefTest(model,H,c,varargin)
%coefTest Linear hypothesis test on coefficients.
%   PVAL = coefTest(GLME) computes the p-value for an F test that all fixed
%   effects parameters in the generalized linear mixed effects model GLME
%   except the intercept are zero.
%
%   PVAL = coefTest(GLME,H) computes the p-value for an F test on the fixed
%   effects part of GLME using the M-by-P matrix, H where P is the number
%   of fixed effects parameters in GLME. Each row of H represents 1
%   contrast, and the columns of H from left to right correspond
%   respectively to the P-by-1 fixed effects vector BETA displayed from top
%   to bottom in the tabular display from the fixedEffects method. The
%   output PVAL is the p-value for an F test that H*BETA = 0. To include
%   contrasts that involve the random effects, use the 'REContrast'
%   name-value pair.
%
%   PVAL = coefTest(GLME,H,C) also specifies a M-by-1 vector C for testing 
%   the hypothesis H*BETA = C.
%
%   PVAL = coefTest(GLME,H,C,PARAM1,VALUE1,...) accepts optional name/value
%   pairs to control the calculation of PVAL.
%
%             Name     Value
%       'DFMethod'     A string that specifies the method to use for
%                      computing the approximate denominator degrees of
%                      freedom DF for the F test. If 'DFMethod' is
%                      'Residual', the DF value is assumed to be equal to
%                      (N-P) where N is the number of observations and P is
%                      the number of fixed effects. If 'DFMethod' is
%                      'None', then the DF value is taken to be infinity.
%                      Default is 'Residual'.
%
%       'REContrast'   A M-by-Q matrix K where Q is the number of random 
%                      effects parameters in GLME. Each row of K represents
%                      1 contrast and the columns of K from left to right
%                      correspond respectively to the Q-by-1 random effects
%                      vector B displayed from top to bottom in the tabular
%                      display from the randomEffects method. The output
%                      PVAL is the p-value for an F test H*BETA + K*B = C.
%
%   [PVAL,F,DF1,DF2] = coefTest(...) also returns the F-statistic F, the
%   numerator degrees of freedom DF1 for F, and the denominator degrees of
%   freedom DF2 for F. DF1 is equal to the number of linearly independent
%   rows in H, or in [H,K] if you use the 'REContrast' name-value pair. The
%   value of DF2 depends on the option selected for 'DFMethod'.
%
%   coefTest uses an approximation to the conditional mean squared error of
%   prediction (CMSEP) of the estimated linear combination of fixed and
%   random effects to compute p-values. This accounts for the uncertainty
%   in the estimation of fixed effects but not the covariance parameters.
%
%   Example: Fit a model with two fixed-effect predictors and a random
%            effect. Test for the significance of the Cylinders term. The
%            p-value is the same as shown in the anova table.
%      load carsmall
%      T = table(MPG,Weight,Model_Year,Cylinders);
%      T.Cylinders = nominal(T.Cylinders);
%      glme = fitglme(T,'MPG ~ Weight + Cylinders + (1|Model_Year)','Distribution','Normal')
%      p = coefTest(glme,[0 0 1 0;0 0 0 1])
%      anova(glme)
%
%      % Test for no difference between 6-cylinder and 8-cylinder cars
%      coefTest(glme,[0 0 1 -1])
%
%   See also anova, coefCI, fixedEffects, randomEffects, covarianceParameters. 

        if nargin < 2
            % Called like: coefTest(model)
            args = {};            
        elseif nargin < 3
            % Called like: coefTest(model,H)
            args = {H};
        elseif nargin < 4
            % Called like: coefTest(model,H,C)
            args = {H,c};
        else
            % Called like: coefTest(model,H,C,PARAM1,VALUE1,...)
            args = [{H,c},varargin];
        end
        [P,F,DF1,DF2] = coefTest@classreg.regr.LinearLikeMixedModel(model,args{:});        

    end % end of coefTest.
    
end

% Other inherited abstract, protected methods.
methods(Access='protected')
    
    function model = fitter(model)% subclass-specific fitting algorithm, must fill in Coefs, CoefficientCovariance, and DFE
                        
        % 1. Save names for the big sparse Zs used to fit a standard GLME.
        model.RandomInfo.ZsColNames = makeSparseZNames(model);
        
        % 2. Fit a GLME model in standard form.        
        model.slme = fitStandardLMEModel(model);
       
        % 3. Fill in DFE and Coefs so that get.NumCoefficients and
        % get.NumEstimatedCoefficients in ParametricRegression work as
        % required. Also fill in CoefficientCovariance.
        model.Coefs = model.slme.betaHat;
        N = model.NumObservations; % only contribution from subset.
        P = length(model.Coefs);
        model.DFE = N - P;
        model.CoefficientCovariance = model.slme.covbetaHat;
        
    end % end of fitter.
    
    % The predict method would normally take a dataset/table array or a
    % matrix containing all variables.  This method exists to allow
    % prediction with a matrix that contains only the required predictor
    % variables without blowing it up to contain all the variables only to
    % then pull it apart to get the design matrix.
    function ypred = predictPredictorMatrix(model,Xpred) %#ok<INUSD>
        ypred = [];
    end % end of predictPredictorMatrix.
    
    function D = get_diagnostics(model,type) %#ok<INUSD>
        D = [];
    end % end of get_diagnostics.

    function L = getlogLikelihood(model)
               
        % L will be the maximized pseudo log likelihood for PL based
        % methods and the maximized log likelihood for ML based methods.
        L = model.ModelCriterion.LogLikelihood;        
        
    end % end of getlogLikelihood.

    function L0 = logLikelihoodNull(model)         %#ok<MANU>

        L0 = NaN;
        
    end % end of logLikelihoodNull.       
    
end

% Other inherited protected methods that we choose to redefine.
methods(Access='protected')
    
    function model = postFit(model)                
                        
        % 1. Set Dispersion, DispersionEstimated and Response.
        model.Dispersion = (model.slme.sigmaHat)^2;
        if ( model.slme.isSigmaFixed == true )
            model.DispersionEstimated = false;
        else
            model.DispersionEstimated = true;
        end
        model.Response = response(model);
        
        % 2. Fill in the log likelihoods.
        model.LogLikelihood     = getlogLikelihood(model);
        model.LogLikelihoodNull = logLikelihoodNull(model);
        
        % 3. Get SSE, SST and SSR.
        [model.SSE,model.SSR,model.SST] = getSumOfSquares(model);
       
        % 4. Get CoefficientNames.
        model.CoefficientNames = getCoefficientNames(model);                
        
        % 5. Save a table of covariance parameters in the object. Computing
        % the CIs on covariance parameters can be expensive - save the
        % covariance parameters table without the CIs.
        [~,~,model.CovarianceTable] = covarianceParameters(model,'WantCIs',false);               
        
        % 6. Fill out BinomSize in model.ObservationInfo for consistency
        % with GeneralizedLinearModel.
        subset = model.ObservationInfo.Subset;
        N = length(subset);
        model.ObservationInfo.BinomSize = ones(N,1);
        model.ObservationInfo.BinomSize(subset) = model.BinomialSize;
        
    end

    function model = selectVariables(model)
        f = model.Formula;
        [~,model.PredLocs] = ismember(f.PredictorNames,f.VariableNames);
        [~,model.RespLoc] = ismember(f.ResponseName,f.VariableNames);
        model = selectVariables@classreg.regr.ParametricRegression(model);
    end
    
end

% Dependent property get methods that we choose to redefine indirectly by
% overloading the methods that will be called in get.property methods in
% superclasses.
methods(Access='protected')
    
    % (1) get.NumCoefficients in ParametricRegression will work if
    % model.Coefs is well defined.
    
    % (2) get.NumEstimatedCoefficients in ParametricRegression will work if
    % model.DFE is well defined.
   
    % (3) Called by get.Coefficients.
    function table = tstats(model)
        [~,~,table] = fixedEffects(model);
    end % end of tstats.
    
    % (4) Called by get.Rsquared. The inherited get_rsquared should work
    % provided model.SSE, model.SST, model.DFE, model.NumObservations,
    % model.LogLikelihood and model.LogLikelihoodNull have sensible values.
    % Not definig a new get_rquared and using the default one.    
    %function rsq = get_rsquared(model)
    %    rsq = [];
    %end % end of get_rsquared.
    
    % (5) Called by get.ModelCriterion.
    function crit = get_modelcriterion(model,type)
                
        crit = modelCriterionLME(model);
        if nargin >= 2
            % type can be one of these: AIC, BIC, LogLikelihood, Deviance
            crit = crit.(type);
        end
        
    end % end of get_modelcriterion.
    
    % (6) Called by get.Residuals.
    function r = get_residuals(model,residualtype)
        
        if nargin < 2
            % (1) Get conditional residuals of various types.
            Raw = residuals(model,'ResidualType','raw');
            Pearson = residuals(model,'ResidualType','pearson');           

            % (2) Assemble into a table.
            r = table(Raw,Pearson,'RowNames',model.ObservationNames);        
        else
            % (1) Ensure that residualtype is a valid 'ResidualType' string.
            residualtype = model.validateResidualType(residualtype);  
            
            % (2) Now get the appropriate conditional residual.
            switch lower(residualtype)                
                case 'raw'
                    r = residuals(model,'ResidualType','raw');
                case 'pearson'
                    r = residuals(model,'ResidualType','pearson');                                 
            end
        end
        
    end % end of get_residuals.
        
    % (7) Called by get.Fitted.
    function yfit = get_fitted(model)     
        
        % (1) Return conditional fitted values - the default in fitted.
        yfit = fitted(model);        
        
    end % end of get_fitted.    
    
end

% No arg constructor.
methods(Access='public', Hidden=true)
   
    function glme = GeneralizedLinearMixedModel(varargin)
%   The GeneralizedLinearMixedModel constructor is not intended to be
%   called directly. Use FITGLME to create a GeneralizedLinearMixedModel by
%   fitting to data.
%
%   See also FITGLME.        
        st = dbstack;
        isokcaller = false;
        if (length(st) >= 2)
            isokcaller = any(strcmpi(st(2).name,{'GeneralizedLinearMixedModel.fit'}));
        end
        if (nargin == 0 && isokcaller == true)
            glme.Formula = classreg.regr.LinearMixedFormula('y ~ -1');
            return;
        end
        error(message('stats:GeneralizedLinearMixedModel:NoConstructor'));
    end
       
end    

% Display related private methods
methods(Access='private')        
    
    function displayHeadLine(model)
        
        % 1. Display headline.
        isLoose = strcmp(get(0,'FormatSpacing'),'loose');
        if (isLoose), fprintf('\n'); end        
        switch lower(model.FitMethod)
            case {'mpl','rempl'}
                headline = [getString(message('stats:GeneralizedLinearMixedModel:Display_headline')),' ','PL'];
            otherwise
                headline = [getString(message('stats:GeneralizedLinearMixedModel:Display_headline')),' ','ML'];
        end
        headline = GeneralizedLinearMixedModel.formatBold(headline);
        fprintf('%s\n',headline);        
        fprintf('\n');
        
    end % end of displayHeadLine.
    
    function displayFormula(model)
                
        % 1. Fit the formula string to the command window width.
        formulaheadline = getString(message('stats:GeneralizedLinearMixedModel:Display_formula'));
        formulaheadline = GeneralizedLinearMixedModel.formatBold(formulaheadline);
        fprintf('%s\n',formulaheadline);
        indent = '    ';
        maxWidth = matlab.desktop.commandwindow.size; 
        maxWidth = maxWidth(1) - 1;
        f = model.Formula;
        fstr = char(f,maxWidth-length(indent));
        disp([indent fstr]);
        fprintf('\n');
        
    end % end of displayFormula.
    
    function displayModelFitStats(model)
        
        % 1. AIC, BIC etc. table.
        modelfitstatsheadline = getString(message('stats:GeneralizedLinearMixedModel:Display_modelfitstats'));
        modelfitstatsheadline = GeneralizedLinearMixedModel.formatBold(modelfitstatsheadline);
        fprintf('%s\n',modelfitstatsheadline);
        crittable = modelCriterionLME(model);        
        crittable = GeneralizedLinearMixedModel.removeTitle(crittable);
        disp(crittable);
        
    end % end of displayModelFitStats.
    
    function displayFixedStats(model)
        
        % 1. Stats table for fixed effects coefficients.
        fixedstatsheadline = getString(message('stats:GeneralizedLinearMixedModel:Display_fixedstats'));
        fixedstatsheadline = GeneralizedLinearMixedModel.formatBold(fixedstatsheadline);
        fprintf('%s\n',fixedstatsheadline);
        ds = model.Coefficients;        
        ds = GeneralizedLinearMixedModel.removeTitle(ds);
        disp(ds);
        
    end % end of displayFixedStats.
    
    function displayCovarianceStats(model)
        
        % 1. Get headline string.
        covariancestatsheadline = getString(message('stats:GeneralizedLinearMixedModel:Display_covariancestats'));
        covariancestatsheadline = GeneralizedLinearMixedModel.formatBold(covariancestatsheadline);
        fprintf('%s\n',covariancestatsheadline);

        % 2. Number of grouping variables
        R = model.GroupingInfo.R;
        
        % 3. If there are random effects, extract the covariance
        % table and display it. Also, display the error covariance.        
        lev = model.GroupingInfo.lev;
        for k = 1:(R+1)
            % 3.1 Name of current group.
            if k > R
                gname = getString(message('stats:GeneralizedLinearMixedModel:String_error'));
                fprintf('%s\n',[getString(message('stats:GeneralizedLinearMixedModel:String_group')),': ',gname]);
            else
                gname = model.GroupingInfo.GNames{k};                
                fprintf('%s\n',[getString(message('stats:GeneralizedLinearMixedModel:String_group')),': ',gname,' (',num2str(lev(k)),' ',getString(message('stats:LinearMixedModel:String_levels')),')']);
            end
            % 3.2 Delete the Group column from covtable{k} if required.
            if isa(model.CovarianceTable{k},'table')
                varnames = ...
                    model.CovarianceTable{k}.Properties.VariableNames;
            elseif isa(model.CovarianceTable{k},'dataset')
                varnames = ...
                    model.CovarianceTable{k}.Properties.VarNames;
            end                
            if any(strcmpi('Group',varnames))
                model.CovarianceTable{k}.Group = [];
            end
            % 3.3 Display covtable{k}.
            ds = model.CovarianceTable{k};
            ds = GeneralizedLinearMixedModel.removeTitle(ds);
            disp(ds);
        end
        
    end % end of displayCovarianceStats.
    
    function displayModelInfo(model)
            
        % 1. Get headling string.
        modelinfoheadline = getString(message('stats:GeneralizedLinearMixedModel:Display_modelinfo'));
        modelinfoheadline = GeneralizedLinearMixedModel.formatBold(modelinfoheadline);
        fprintf('%s\n',modelinfoheadline);
        
        % 2. Number of observations.
        N = model.slme.N;
        
        % 3. Number of fixed effects.
        p = model.slme.p;
        
        % 4. Number of random effects.
        q = model.slme.q;        

        % 5. Number of covariance parameters.
        if model.slme.isSigmaFixed
            ncov = model.slme.Psi.NumParametersExcludingSigma;
        else
            ncov = model.slme.Psi.NumParametersExcludingSigma + 1;
        end
        
        % 6. Distribution.
        distribution = model.convertFirstCharToUpper(model.Distribution);        
        
        % 7. Link.
        linkname = model.convertFirstCharToUpper(model.Link.Name);
        
        % 8. FitMethod.
        switch lower(model.FitMethod)
            case 'mpl'
                fitmethod = 'MPL';
            case 'rempl'
                fitmethod = 'REMPL';
            case 'approximatelaplace'
                fitmethod = 'ApproximateLaplace';
            case 'laplace'
                fitmethod = 'Laplace';
        end
        
        % 9. Show the N, p, q and ncov in a dataset.
        indent = '    ';
        fprintf('%-35s %6d\n',[indent,getString(message('stats:GeneralizedLinearMixedModel:ModelInfo_numobs'))],N);
        fprintf('%-35s %6d\n',[indent,getString(message('stats:GeneralizedLinearMixedModel:ModelInfo_fecoef'))],p);
        fprintf('%-35s %6d\n',[indent,getString(message('stats:GeneralizedLinearMixedModel:ModelInfo_recoef'))],q);
        fprintf('%-35s %6d\n',[indent,getString(message('stats:GeneralizedLinearMixedModel:ModelInfo_covpar'))],ncov);
        fprintf('%-35s %-6s\n',[indent,getString(message('stats:GeneralizedLinearMixedModel:ModelInfo_distribution'))],distribution);
        fprintf('%-35s %-6s\n',[indent,getString(message('stats:GeneralizedLinearMixedModel:ModelInfo_link'))],linkname);
        fprintf('%-35s %-6s\n',[indent,getString(message('stats:GeneralizedLinearMixedModel:ModelInfo_fitmethod'))],fitmethod);
        fprintf('\n');        
        
    end % end of displayModelInfo.
    
end

% Private utility methods.
methods (Access={?classreg.regr.LinearLikeMixedModel})
                                  
    function crittable = modelCriterionLME(model)
%modelCriterionLME - Compute table containing model criterion info.
%   crittable = modelCriterionLME(model) takes GeneralizedLinearMixedModel
%   object model and computes a table of model criterion such as AIC, BIC,
%   LogLikelihood and Deviance.         
    
        % 1. Call modelCriterion method in the implementation class.
        crit = modelCriterion(model.slme);    
        
        % 2. Put crit into a table with desired variable names.
        crittable = table(crit.AIC, crit.BIC, crit.logLik, crit.Deviance,...
            'VariableNames',{'AIC' 'BIC' 'LogLikelihood' 'Deviance'});        
        
        % 3. Put a title on the output table.
        ttl = getString(message('stats:LinearMixedModel:Title_modelfitstats'));
        crittable = classreg.regr.lmeutils.titleddataset(crittable,ttl);                
        
    end % end of modelCriterionLME.
        
    function w = getCombinedWeights(model,reduce)
%getCombinedWeights - Get the effective weights for this model.        
%   w = getCombinedWeights(model,reduce) takes an object model of type
%   GeneralizedLinearMixedModel and returns w, the effective observation
%   weights vector. If reduce is true then the length of w equals the
%   number of observations actually used in the fit. If reduce is false
%   then the length of w equals the total number of observations in the
%   input data even if not used in the fit.

        % (1) Get current observation weights.
        w = model.ObservationInfo.Weights;
        
        % (2) If required, return weights only for the observations
        % actually used in the fit.
        if nargin<2 || reduce
            subset = model.ObservationInfo.Subset;
            w = w(subset);
        end
        
    end % end of getCombinedWeights.
    
    function slme = fitStandardLMEModel(model)
%fitStandardLMEModel - Fits a GLME model in standard form.
%   slme = fitStandardLMEModel(model) takes an object model of type
%   GeneralizedLinearMixedModel and returns a fitted model slme of class
%   StandardGeneralizedLinearMixedModel that is ready for stats.

        % 1. Make big sparse Zs matrix and covariance matrix object Psi.
        Z = model.RandomInfo.Z;
        q = model.RandomInfo.q;
        lev = model.GroupingInfo.lev;
        Gid = model.GroupingInfo.Gid;
        N = model.NumObservations;
        Zs = GeneralizedLinearMixedModel.makeSparseZ(Z,q,lev,Gid,N);
        Psi = makeCovarianceMatrix(model);

        % 2. Exclude observations, NaNs or Infs. Already done.

        % 3. Get reduced weights and fixed effects design matrix.
        reduce = true;
        w = getCombinedWeights(model,reduce);        
        X = model.FixedInfo.X;                     

        % 4. Fit a StandardGeneralizedLinearMixedModel and make it ready for stats.
        dofit = true;
        dostats = true;        
        args = {'Distribution',model.Distribution,...
            'BinomialSize',model.BinomialSize,...
            'Link',model.Link,...
            'Offset',model.Offset,...
            'DispersionFlag',model.DispersionFlag,...
            'Weights',w,...
            'Optimizer',model.Optimizer,...
            'OptimizerOptions',model.OptimizerOptions,...
            'InitializationMethod',model.StartMethod,...
            'CheckHessian',model.CheckHessian,...
            'PLIterations',model.PLIterations,...
            'PLTolerance',model.PLTolerance,...
            'MuStart',model.MuStart,...
            'InitPLIterations',model.InitPLIterations,...
            'EBMethod',model.EBMethod,...
            'EBOptions',model.EBOptions,...
            'CovarianceMethod',model.CovarianceMethod,...
            'UseSequentialFitting',model.UseSequentialFitting,...
            'ShowPLOptimizerDisplay',model.ShowPLOptimizerDisplay};
        slme = classreg.regr.lmeutils.StandardGeneralizedLinearMixedModel(X,model.y,Zs,Psi,model.FitMethod,dofit,dostats,args{:});

    end % end of fitStandardLMEModel.
            
    function [SSE,SSR,SST] = getSumOfSquares(model)
%getSumOfSquares - Computes SSE, SSR and SST.
%   [SSE,SSR,SST] = getSumOfSquares(model) takes an object model of type
%   GeneralizedLinearMixedModel and computes SSE, SSR and SST using
%   conditional fitted values. Only observations actually used in the fit
%   should be used for calculating SSE, SSR and SST.
%
%   Suppose w is a vector of observation weights. w is the product of prior
%   weights and binomial size for binomial distribution and is equal to
%   prior weights for other distributions. If Y is the response vector and
%   F is a vector of conditional fitted values then:
%
%   SSE = sum( weff.* (Y-F).^2 );
%   SSR = sum( weff.* (F-meanw(F)).^2 );
%   SST = SSE + SSR;
%  
%   where weff(i) = w(i)/v(i) and v is the variance function evaluated at
%   F. meanw(F) is the weighted mean of F using weights weff.

        % 1. Subset of observations actually used in the fit.
        subset = model.ObservationInfo.Subset;
        
        % 2. Get *conditional* fitted values.
        F = fitted(model);
        F = F(subset);
        
        % 3. Get response vector.
        Y = response(model);
        Y = Y(subset);
        
        % 4. Get the weight vector.
        w = model.ObservationInfo.Weights;
        w = w(subset);
        if strcmpi(model.Distribution,'binomial')
            w = w .* model.BinomialSize;
        end
        
        % 5. Compute "effective" weights.
        vfun = model.slme.VarianceFunction.VarianceFunction;
        v = vfun(F);        
        weff = w./v;                
        
        % 6. Weighted mean of F.
        F_mean_w = sum(weff.*F)/sum(weff);
        
        % 7. Compute SSE.
        SSE = sum(weff.*((Y - F).^2));   
        
        % 8. Compute SSR.
        SSR = sum(weff.*((F - F_mean_w).^2));
        
        % 9. Compute SST.
        SST = SSE + SSR;                                   

    end % end of getSumOfSquares.
    
    function coefnames = getCoefficientNames(model)
%getCoefficientNames - Get names associated with the rows in Coefficients.
%   coefnames = getCoefficientNames(model) returns a cell array of string
%   containing a label for coefficients that appear in rows of Coefficients
%   property.

        % (1) model.FixedInfo.XColNames is what you want.
        coefnames = model.FixedInfo.XColNames;

    end % end of getCoefficientNames.    
           
    function np = getTotalNumberOfParameters(model)
%getTotalNumberOfParameters - Get the total number of parameters in model.        
%   np = getTotalNumberOfParameters(model) returns the total number of 
%   fixed effects parameters plus the total number of covariance parameters 
%   plus 1 (if the dispersion parameter is not fixed).
       
        % 1. Get number of fixed effects parameters.
        numfixpar = model.slme.p;
        
        % 2. Get number of covariance parameters excluding dispersion.
        numcovpar = model.slme.Psi.NumParametersExcludingSigma;
        
        % 3. Get np.
        if model.slme.isSigmaFixed 
            np = numfixpar + numcovpar;
        else
            np = numfixpar + numcovpar + 1;
        end
        
    end % end of getTotalNumberOfParameters.
    
    function w = getEffectiveObservationWeights(model,reduce)
%w = getEffectiveObservationWeights(model,reduce) combines the prior 
%   weights with binomial size if required and gets the effective
%   observation weights. For most distributions the prior weights are the
%   same as effective weights. The one exception is the binomial
%   distribution for which the prior weights need to be multiplied by the
%   binomial size. If reduce is true, only effective observation weights
%   corresponding to the subset of observations actually used in the fit is
%   returned.
            
            % 1. Get the prior weights.
            w = model.ObservationInfo.Weights;
            
            % 2. Get the binomial size.
            binomsize = model.ObservationInfo.BinomSize;

            % 3. Create "effective" observation weights.
            if strcmpi(model.Distribution,'binomial')
                w = w.*binomsize;
            end
            
            % 4. Return either reduced or full w.
            if nargin<2 || reduce
                subset = model.ObservationInfo.Subset;
                w = w(subset);
            end
            
    end % end of getEffectiveObservationWeights.
    
end

% Private static utility methods.
methods(Static,Access='private')
                                             
    function strout = convertFirstCharToUpper(strin)
%strout = convertFirstCharToUpper(strin) takes a length M string strin and
%   converts its first character to upper case. strout is the modified
%   string.

        % 1. Do nothing if strin is empty.
        if isempty(strin)
            strout = strin;
            return;
        end
        
        % 2. Make first element of strin uppercase.
        strout    = strin;
        strout(1) = upper(strout(1));
        
    end % end of convertFirstCharToUpper.
    
end

methods(Static,Access='protected')
    
    function checkNestingRequirement(smallModel,bigModel,smallModelName,bigModelName)
%checkNestingRequirement - Checks if smallModel is nested in bigModel.       
%  checkNestingRequirement(smallModel,bigModel,smallModelName,bigModelName)
%  takes two GeneralizedLinearMixedModel objects, smallModel and bigModel.
%  smallModel is the potentially smaller model that is nested in the bigger
%  model bigModel. smallModel has name smallModelName and bigModel has name
%  bigModelName. We try to check if smallModel is really nested in bigModel
%  or not. This technique is not fool proof since we do not check the
%  covariance structure of random effects but it is better than no check.
%
%   What is checked?
%  
%   (1) smallModel and bigModel have been fit using the same FitMethod. If
%   the models have been fit using 'rempl', they cannot be compared.
%
%   (2) smallModel and bigModel are being modelled using the same
%   distribution.
% 
%   (3) smallModel and bigModel use the same link function.
%
%   (4) smallModel and bigModel have the same effective observation weights 
%   vector (used in fit).
%
%   (5) smallModel and bigModel have the same response vector used in fit.
%
%   (6) The maximized log likelihood of bigModel must be >= the maximized 
%   log likelihood of the smallModel.
%
%   Suppose Xsmall is the fixed effects design matrix of smallModel in
%   standard form and Xbig is the fixed effects design matrix of bigModel
%   in standard form actually used in the fit. 
%
%   (7) The span of Xbig must contain Xsmall.
%
%   Suppose Zsmall is the overall random effects design matrix of
%   smallModel in standard form and Zbig is the overall random effects
%   design matrix of bigModel in standard form actually used in the fit.
%
%   (8) The span of Zbig must contain Zsmall.

        % 1. smallModel/bigModel are GeneralizedLinearMixedModel objects?
        assert( isa(smallModel,'GeneralizedLinearMixedModel') );
        assert( isa(  bigModel,'GeneralizedLinearMixedModel') );
        
        % 2. Ensure that smallModelName and bigModelName are strings.
        assert( internal.stats.isString(smallModelName) );
        assert( internal.stats.isString(  bigModelName) );       
                
        % 3. Do the nesting checks.
        
        % 3.1 Have smallModel and bigModel been fit using the same FitMethod
        % that is not equal to 'REMPL'.
        fitmethodsmall = smallModel.FitMethod;
        fitmethodbig   =   bigModel.FitMethod;
        isok = strcmpi(fitmethodsmall,fitmethodbig) & (strcmpi(fitmethodsmall,'rempl') == false);
        assertThat(isok,'stats:GeneralizedLinearMixedModel:NestingCheck_fitmethod',smallModelName,bigModelName);        
        
        % 3.2 Ensure that smallModel and bigModel are using the same
        % distribution for the response.
        distrsmall = smallModel.Distribution;
        distrbig   =   bigModel.Distribution;
        assertThat(strcmpi(distrsmall,distrbig),'stats:GeneralizedLinearMixedModel:NestingCheck_distribution',smallModelName,bigModelName);        
                        
        % 3.3 smallModel and bigModel must use the same link function.
        linksmall = smallModel.Link.Name;
        linkbig   =   bigModel.Link.Name;
        assertThat(isequal(linksmall,linkbig),'stats:GeneralizedLinearMixedModel:NestingCheck_link',smallModelName,bigModelName);
        
        % 3.4 For 'MPL', if both models have been fitted using the normal
        % distribution and identity link then it is OK to compare them.
        % Otherwise, warn that we are using the maximized log likelihood of
        % pseudo data from final PL iteration.
        if strcmpi(fitmethodsmall,'mpl')
            oktocompare = any(strcmpi(distrsmall,{'normal','gaussian'})) && strcmpi(linksmall,'identity');
            if ~oktocompare
                warning(message('stats:GeneralizedLinearMixedModel:LRTUsingPseudoData'));
            end
        end
        
        % 3.5 smallModel and bigModel must have the same "effective"
        % observation weights.
        reduce = true;
        wsmall = getEffectiveObservationWeights(smallModel,reduce);
        wbig   = getEffectiveObservationWeights(  bigModel,reduce);
        assertThat(max(abs(wsmall-wbig)) <= sqrt(eps),'stats:GeneralizedLinearMixedModel:NestingCheck_weights',smallModelName,bigModelName);
        
        % 3.6 Ensure that smallModel and bigModel have the same response
        % vector, actually used in the fit.
        ysmall = smallModel.y;
        ybig   =   bigModel.y;  
        assertThat(isequaln(ysmall,ybig),'stats:GeneralizedLinearMixedModel:NestingCheck_response',smallModelName,bigModelName);
        
        % 3.7  Maximized log likelihood of bigModel >= that of smallModel?
        logliksmall = smallModel.LogLikelihood;
        loglikbig   =   bigModel.LogLikelihood;
        assertThat(loglikbig >= logliksmall,'stats:GeneralizedLinearMixedModel:NestingCheck_loglik',bigModelName,smallModelName);
        
        % 3.8 The span of Xbig must contain Xsmall.
        Xsmall = smallModel.FixedInfo.X;
        Xbig   =   bigModel.FixedInfo.X;
        assertThat(GeneralizedLinearMixedModel.isMatrixNested(Xsmall,Xbig),'stats:GeneralizedLinearMixedModel:NestingCheck_nestedspanX',smallModelName,bigModelName);
                
        % 3.9 The span of Zbig must contain Zsmall.                                                           
        Zsmall = smallModel.slme.Z;
        Zbig   =   bigModel.slme.Z;
        assertThat(GeneralizedLinearMixedModel.isMatrixNested(Zsmall,Zbig),'stats:GeneralizedLinearMixedModel:NestingCheck_nestedspanZ',smallModelName,bigModelName);
        
    end % end of checkNestingRequirement.
    
end

% Validation methods that we are required to define by subclassing
% LinearLikeMixedModel.
methods(Static,Access='protected')
    
    function fitmethod = validateFitMethod(fitmethod)
%validateFitMethod - Validate the 'FitMethod' parameter.
%   fitmethod = validateFitMethod(fitmethod) takes a potential 'FitMethod'
%   parameter and either returns the validated parameter fitmethod or 
%   throws an error message. 
%
%   What is checked?
%
%   (1) fitmethod is a string defined in AllowedFitMethods property.
%
%   If (1) is *not* satisfied, then an error message is thrown.
        
        fitmethod = internal.stats.getParamVal(fitmethod,...
            GeneralizedLinearMixedModel.AllowedFitMethods,'FitMethod');

    end % end of validateFitMethod.

    function w = validateWeights(w,N,distribution)
%validateWeights - Validates the 'Weights' parameter.        
%   w = validateWeights(w,N,distribution) takes a potential 'Weights'
%   parameter w, an integer N specifying the expected size of w and a
%   string distribution specifying the conditional distribution of y given
%   b. Either the output w is a validated column vector of weights or an
%   error is thrown.
%
%   What is checked?
%
%   (1) w is a numeric, real, vector of length N.
%
%   (2) All elements of w are not-NaN, >= 0 and < Inf.
%
%   (3) If distribution is 'binomial' or 'poisson' then all elements of w
%       must be positive integers.
%
%   If any of (1), (2) or (3) do not hold, then an error message is thrown. 

        % 1. N must be >= 0 and it must be an integer.
        assert(N >= 0 & internal.stats.isScalarInt(N));        
        
        % 2. (a) w must be a numeric, real, vector with length(w) == N. 
        %    (b) All elements of w are non-NaN, >= 0 and < Inf.
        strN = num2str(N);
        assertThat(isnumeric(w) & isreal(w) & isvector(w) & length(w)==N,'stats:GeneralizedLinearMixedModel:BadWeights',strN,strN);
        assertThat(      all( ~isnan(w) & w >= 0 & w < Inf )            ,'stats:GeneralizedLinearMixedModel:BadWeights',strN,strN);        
   
        % 3. For 'binomial' and 'poisson' distributions, elements of w must
        % be positive integers.
        if any(strcmpi(distribution,{'binomial','poisson'}))
            assertThat(internal.stats.isIntegerVals(w,1),'stats:GeneralizedLinearMixedModel:BadWeights',strN,strN);
        end
        
        % 4. Make w into a column vector if required.
        if size(w,1) == 1
            w = w';
        end
        
    end % end of validateWeights.

    function [optimizer,optimizeroptions] = ...
            validateOptimizerAndOptions(optimizer,optimizeroptions)
%validateOptimizerAndOptions - Validate 'Optimizer'/'OptimizerOptions'.
%   [optimizer,optimizeroptions] = validateOptimizerAndOptions(optimizer,optimizeroptions)
%   takes a potential string optimizer and a structure or an object
%   optimizeroptions for the optimizer and returns the validated values.
%
%   What is checked?
%
%   (1) optimizer is a string that occurs in the cell array AllowedOptimizers.
%
%   (2) If optimizer is 'quasinewton' or 'fminsearch', optimzeroptions can 
%   either be empty or a struct.
%
%   (3) If optimizer is 'fminunc', optimizeroptions can either be empty or
%   an object of class optim.options.SolverOptions.
%
%   (4) When using 'fminunc', check for Optimization Toolbox license.
%
%   (5) If optimizeroptions is empty []:
%           (a) If optimizer is 'quasinewton' then optimizeroptions =
%           statset('linearmixedmodel').
%           (b) If optimizer is 'fminunc' then optimizeroptions =
%           optimoptions('fminunc').
%           (c) If optimizer is 'fminsearch' then optimizeroptions = 
%           optimset('fminsearch').

        % 1. Ensure that optimizer is sensible.
        optimizer = internal.stats.getParamVal(optimizer,...
            GeneralizedLinearMixedModel.AllowedOptimizers,'Optimizer');        
        
        % 2. Check for license if required.
        switch lower(optimizer)
            case {'quasinewton','fminsearch'}                                                            
                % 2.1 No license check needed.                
            case {'fminunc'}                
                % 2.2 Check for license when using fminunc.
                if ( license('test','optimization_toolbox') == false )
                    error(message('stats:GeneralizedLinearMixedModel:LicenseCheck_fminunc'));                    
                end
        end
        
        % 3. Do the following:
        %   (a) Ensure that optimizeroptions is either empty or a struct or
        %   object created using statset/optimset or optimoptions depending
        %   on the string optimizer and
        %
        %   (b) Create a completely filled out optimizeroptions structure 
        %   or object as required. Either use the defaults or merge the 
        %   supplied options with the defaults.
        switch lower(optimizer)
            case 'quasinewton'                                                                            
                % 3.1 optimizeroptions can be empty or a struct, otherwise 
                % we have an error.
                assertThat(isempty(optimizeroptions) || isstruct(optimizeroptions),'stats:GeneralizedLinearMixedModel:OptimizerOptions_qn'); 
                
                % 3.2 Filled optimizeroptions.
                dflts = statset('linearmixedmodel');
                if isempty(optimizeroptions)
                    optimizeroptions = dflts;
                else
                    optimizeroptions = statset(dflts,optimizeroptions);
                end
                
            case {'fminunc'}                               
                % 3.1 optimizeroptions can be empty or of class
                % optim.options.SolverOptions, otherwise we have an error.
                assertThat(isempty(optimizeroptions) || isa(optimizeroptions,'optim.options.SolverOptions'),'stats:GeneralizedLinearMixedModel:OptimizerOptions_fminunc');        
                
                % 3.2 Filled optimizeroptions.
                if isempty(optimizeroptions)
                    optimizeroptions = optimoptions('fminunc');
                    optimizeroptions.Algorithm = 'quasi-newton';
                else
                    optimizeroptions = optimoptions('fminunc',optimizeroptions);
                end
                
            case 'fminsearch'
                % 3.1 optimizeroptions can be empty or a struct, otherwise 
                % we have an error.
                assertThat(isempty(optimizeroptions) || isstruct(optimizeroptions),'stats:GeneralizedLinearMixedModel:OptimizerOptions_fminsearch'); 
                
                % 3.2 Filled optimizeroptions.
                dflts = optimset('fminsearch');
                if isempty(optimizeroptions)
                    optimizeroptions = dflts;
                else
                    optimizeroptions = optimset(dflts,optimizeroptions);
                end
        end
        
    end % end of validateOptimizerAndOptions.    
    
    function startmethod = validateStartMethod(startmethod)
%validateStartMethod - Validates the 'StartMethod' parameter.
%   startmethod = validateStartMethod(startmethod) accepts a potential
%   value of 'StartMethod' parameter and validates it. If not valid, an
%   error message is thrown.
%
%   What is checked?
%
%   (1) startmethod is a string that occurs in the cell array AllowedStartMethods.
%
%   If (1) is *not* satisfied, then an error message is thrown.

        startmethod = ...
            internal.stats.getParamVal(startmethod,GeneralizedLinearMixedModel.AllowedStartMethods,'StartMethod');     
            
    end % end of validateStartMethod.

    function dfmethod = validateDFMethod(dfmethod)
%validateDFMethod - Validates the degrees of freedom method.        
%   dfmethod = validateDFMethod(dfmethod) accepts a potential degrees of 
%   freedom and validates it. If not valid, an error message is thrown.
%
%   What is checked?
%
%   (1) dfmethod is a string contained in AllowedDFMethods.
%
%   If (1) is not true, an error message is thrown.

        dfmethod = internal.stats.getParamVal(dfmethod,...
            GeneralizedLinearMixedModel.AllowedDFMethods,'DFMethod');     
        
    end % end of validateDFMethod.

    function residualtype = validateResidualType(residualtype)
%validateResidualType - Validates the 'ResidualType' parameter.
%   residualtype = validateResidualType(residualtype) takes a potential
%   value of 'ResidualType' parameter and validates it. If not valid, an
%   error is thrown.
%
%   What is checked?
%
%   (1) residualtype must be a string contained in AllowedResidualTypes.
%
%   If (1) is not satisfied, an error message is thrown.

        residualtype = internal.stats.getParamVal(residualtype,...
            GeneralizedLinearMixedModel.AllowedResidualTypes,'ResidualType');
        
    end % end of validateResidualType.                                                             

    function verbose = validateVerbose(verbose)
%validateVerbose - Validate the parameter 'Verbose'.
%   verbose = validateVerbose(verbose) takes a potential value verbose for
%   the 'Verbose' input and validates it. If not valid, an error message is
%   thrown.
%
%   What is checked?
%
%   (1) verbose must be a scalar logical (true or false). 
% 
%   or
%
%   (2) verbose must be an integer - either 0, 1 or 2.
%
%   If neither (1) nor (2) holds, then an error message is thrown.
        
        % 1. Do we have a valid scalar logical?
        isvalidscalarlogical = isscalar(verbose) && islogical(verbose);
        
        % 2. Do we have a valid scalar int?
        isvalidscalarint = internal.stats.isScalarInt(verbose,0,2);
        
        % 3. Is verbose valid?
        isok = isvalidscalarlogical || isvalidscalarint;
        if ~isok
            error(message('stats:GeneralizedLinearMixedModel:BadVerbose'));
        end
        
        % 4. If verbose is logical, convert it into double.
        if isvalidscalarlogical
            verbose = double(verbose);
        end
        
        % 5. At this point, verbose is either 0, 1 or 2.        

    end % end of validateVerbose.

    function checkhessian = validateCheckHessian(checkhessian)
%validateCheckHessian - Validate the parameter 'CheckHessian'.
%   checkhessian = validateCheckHessian(checkhessian) takes a potential
%   value checkhessian for the 'CheckHessian' flag and validates it. If not
%   valid, an error message is thrown.
%
%   What is checked?
%
%   (1) checkhessian must be a scalar logical (true or false).
%
%   If (1) is not satisfied, then an error message is thrown.
        
        checkhessian = ...
            GeneralizedLinearMixedModel.validateLogicalScalar(checkhessian,'stats:GeneralizedLinearMixedModel:BadCheckHessian');

    end % end of validateCheckHessian.    
    
end

% Protected, static utility methods for input validation.
methods(Static,Access='protected')    
    
    function distribution = validateDistribution(distribution)
%validateDistribution - Validate the 'Distribution' parameter.
%   distribution = validateDistribution(distribution) takes a potential
%   value for 'Distribution' and returns the validated value or throws an
%   error message.
%
%   What is checked?
%
%   (1) distribution must be a string defined in AllowedDistributions
%   property.
%
%   If (1) is *not* satisfied, then an error message is thrown.

        distribution = internal.stats.getParamVal(distribution,...
            GeneralizedLinearMixedModel.AllowedDistributions,'Distribution');
        
    end % end of validateDistribution.
                                                
    function binomialsize = validateBinomialSize(binomialsize,N)
%binomialsize = validateBinomialSize(binomialsize,N) takes a potential 
%   value for the 'BinomialSize' parameter and an integer N and returns the 
%   validated N-by-1 vector binomialsize.
%
%   What is checked?
%
%   (1) binomialsize must be a scalar or a vector of length N.
%
%   (2) If scalar, binomialsize must be an integer > 0. If vector, all
%   elements of binomialsize must be integers > 0.
           
            % 1. Ensure that N is a scalar integer.
            assert(internal.stats.isScalarInt(N));
            
            % 2. Ensure that binomialsize is a N-by-1 vector of positive
            % integers.
            isintvals = internal.stats.isIntegerVals(binomialsize,1);            
            isokscalar = isintvals & isscalar(binomialsize);            
            isokvector = isintvals & isvector(binomialsize) & length(binomialsize) == N;                        
            assertThat(isokscalar | isokvector,'stats:GeneralizedLinearMixedModel:BadBinomialSize',num2str(N));
            
            % 3. If scalar, make binomialsize into a N-by-1 vector.
            if isokscalar
                binomialsize = binomialsize*ones(N,1);
            end
            
            % 4. If binomialsize is a row vector, make it into a column
            % vector.
            if size(binomialsize,1) == 1
                binomialsize = binomialsize';
            end
            
    end % end of validateBinomialSize.           
    
    function linkSpec = defaultLink(distribution)
%linkSpec = defaultLink(distribution) takes a string distribution that 
%   represents one of the allowed distributions and returns the default
%   (canonical) link for the distribution.                 

            switch lower(distribution)            
                case {'normal','gaussian'}
                    linkSpec = 'identity';
                case 'binomial'
                    linkSpec = 'logit';
                case 'poisson'
                    linkSpec = 'log';
                case 'gamma'
                    linkSpec = -1;
                case {'inverse gaussian','inversegaussian'}
                    linkSpec = -2;
            end
                
    end % end of defaultLink.

    function offset = validateOffset(offset,N)
%offset = validateOffset(offset,N) takes a vector offset and an integer N, 
%   and returns the validated offset vector. offset must be a numeric, real
%   vector of size N-by-1.

            % 1. Ensure that N is a scalar integer.
            assert(internal.stats.isScalarInt(N));
            
            % 2. Ensure that offset is a N-by-1 numeric, real vector.
            isok = isnumeric(offset) & isreal(offset) & ...
                isvector(offset) & length(offset) == N;            
            assertThat(isok,'stats:GeneralizedLinearMixedModel:BadOffset',num2str(N));
            
            % 3. If offset is a row vector, make it into a column vector.
            if size(offset,1) == 1
                offset = offset';
            end            
            
    end % end of validateOffset.    

    function dispersionflag = validateDispersionFlag(dispersionflag)
%validateDispersionFlag - Validate the parameter 'DispersionFlag'.
%   dispersionflag = validateDispersionFlag(dispersionflag) takes a
%   potential value dispersionflag for 'DispersionFlag' and validates it.
%   If not valid, an error message is thrown.
%
%   What is checked?
%
%   (1) dispersionflag must be a scalar logical (true or false).
%
%   If (1) is not satisfied, then an error message is thrown.
        
        dispersionflag = ...
            GeneralizedLinearMixedModel.validateLogicalScalar(dispersionflag,'stats:GeneralizedLinearMixedModel:BadDispersionFlag');

    end % end of validateDispersionFlag.    
        
    function pliterations = validatePLIterations(pliterations)
%pliterations = validatePLIterations(pliterations) takes a potential value 
%   for the number of PL iterations and returns the validated value.
%
%   What is checked?
%
%   (1) pliterations must be a positive scalar integer.
            
            isok = isscalar(pliterations) & ...
                internal.stats.isIntegerVals(pliterations,1);            
            assertThat(isok,'stats:GeneralizedLinearMixedModel:BadPLIterations');            
            
    end % end of validatePLIterations.
    
    function pltolerance = validatePLTolerance(pltolerance)
%pltolerance = validatePLTolerance(pltolerance) takes a potential value for 
%   the tolerance for PL iterations and returns the validated value.
%
%   What is checked?
%
%   (1) pltolerance must be a numeric, real scalar.
            
            isok = isscalar(pltolerance) & ...
                isnumeric(pltolerance) & isreal(pltolerance);
            assertThat(isok,'stats:GeneralizedLinearMixedModel:BadPLTolerance');           
            
    end % end of validatePLTolerance.
        
    function mustart = validateMuStart(mustart,distribution,N)
%mustart = validateMuStart(mustart,distribution,N) takes a potential value 
%   of sglme.MuStart and validates it. Legal values of elements in mustart 
%   depend on the value of distribution as follows:
%
%                 Distribution          Legal values
%                 'normal'              (-Inf,Inf)
%                 'binomial'            (0,1)
%                 'poisson'             (0,Inf)
%                 'gamma'               (0,Inf)
%                 'inverse gaussian'    (0,Inf)      
%
%   mustart is required to be a numeric, real vector of size N-by-1. Empty
%   values of mustart are acceptable.
                        
            if ~isempty(mustart)
                % 1. Non-empty mustart must be N-by-1 numeric, real vector.
                sizeok = isnumeric(mustart) & isreal(mustart) & isvector(mustart) & (length(mustart) == N);
                
                % 2. Impose distribution specific constraints.
                switch lower(distribution)
                    case {'binomial'}
                        isok = sizeok & all(mustart > 0 & mustart < 1);                        
                    case {'poisson'}
                        isok = sizeok & all(mustart > 0 & mustart < Inf);
                    case {'gamma'}
                        isok = sizeok & all(mustart > 0 & mustart < Inf);
                    case {'inverse gaussian','inversegaussian'}
                        isok = sizeok & all(mustart > 0 & mustart < Inf);
                    case {'normal','gaussian'}
                        isok = sizeok & all(mustart > -Inf & mustart < Inf);                        
                end
                assertThat(isok,'stats:GeneralizedLinearMixedModel:BadMuStart',num2str(N));                
                
                % 3. If mustart is a row vector, make it into a column
                % vector.
                if size(mustart,1) == 1
                    mustart = mustart';
                end
            end                
            
    end % end of validateMuStart.
    
    function initpliterations = validateInitPLIterations(initpliterations)
%initpliterations = validateInitPLIterations(initpliterations) takes a 
%   potential value of initial number of PL iterations and validates it.
%
%   What is checked?
%
%   (1) initpliterations must be a scalar positive integer.

            isok = isscalar(initpliterations) & ...
                internal.stats.isIntegerVals(initpliterations,1);
            assertThat(isok,'stats:GeneralizedLinearMixedModel:BadInitPLIterations');            

    end % end of validateInitPLIterations.    
    
    function [ebmethod,eboptions] = validateEBParameters(ebmethod,eboptions)
%[ebmethod,eboptions] = validateEBParameters(ebmethod,eboptions) takes
%   ebmethod - a potential value for 'EBMethod' and eboptions - a potential
%   value for 'EBOptions' and returns the validated values.
%
%   What is checked?
%
%   (1) ebmethod must be defined in AllowedEBMethods. If ebmethod is
%   'auto', it is changed to 'default'.
%
%   (2) If eboptions is [] and ebmethod is 'fsolve' then eboptions is set 
%   to optimoptions('fsolve') with default values for TolFun, TolX, MaxIter
%   and Display taken from dfltEBOptions (see below).
%
%   (3) If eboptions is [] and ebmethod is not 'fsolve' then eboptions is
%   set to dfltEBOptions (see below).
%
%   (4) If eboptions is not [] and ebmethod is 'fsolve' then eboptions must
%   be of class optimoptions('fsolve').
%
%   (5) If eboptions is not [] and ebmethod is not 'fsolve' then eboptions
%   must be a structure. eboptions is combined with dfltEBOptions to fill
%   out unspecified values.

            % 1. Create default options.
            dfltEBOptions = statset('TolFun',1e-6,'TolX',1e-8,'MaxIter',100,'Display','off');                                    
            
            % 2. Validate ebmethod. Change ebmethod 'auto' to 'default'.
            ebmethod = internal.stats.getParamVal(ebmethod,GeneralizedLinearMixedModel.AllowedEBMethods,'EBMethod');
            if strcmpi(ebmethod,'auto')
                ebmethod = 'default';
            end
            
            % 3. Validate eboptions.
            switch lower(ebmethod)                
                case 'fsolve'                    
                    % 3.1 Check for license when using fsolve.
                    if ( license('test','optimization_toolbox') == false )
                        error(message('stats:GeneralizedLinearMixedModel:LicenseCheck_fsolve'));
                    end
                    
                    if isempty(eboptions)
                        % 3.2 Create the right object and copy over values
                        % from dfltEBOptions.
                        eboptions = optimoptions('fsolve');
                        eboptions.TolFun  = dfltEBOptions.TolFun;
                        eboptions.TolX    = dfltEBOptions.TolX;
                        eboptions.MaxIter = dfltEBOptions.MaxIter;
                        eboptions.Display = dfltEBOptions.Display;
                    else
                        % 3.3 We must have the right object.
                        assertThat(isa(eboptions,'optim.options.Fsolve'),'stats:GeneralizedLinearMixedModel:BadEBOptions');
                    end
                otherwise
                    % not fsolve
                    if isempty(eboptions)
                        % 3.3 Just return the default options.
                        eboptions = dfltEBOptions;
                    else
                        % 3.4 eboptions must be a structure.
                        assertThat(isstruct(eboptions),'stats:GeneralizedLinearMixedModel:BadEBOptions');
                        % 3.5 Combine eboptions with dfltEBOptions.
                        eboptions = statset(dfltEBOptions,eboptions);
                    end
            end
            
    end % end of validateEBParameters.
        
    function covariancemethod = validateCovarianceMethod(covariancemethod)
%covariancemethod = validateCovarianceMethod(covariancemethod) takes a 
%   potential value for 'CovarianceMethod' parameter and returns the
%   validated value.
%
%   What is checked?
%
%   (1) covariancemethod must be defined in AllowedCovarianceMethods.

        covariancemethod = internal.stats.getParamVal(covariancemethod,...
            GeneralizedLinearMixedModel.AllowedCovarianceMethods,'CovarianceMethod'); 

    end % end of validateCovarianceMethod.
    
    function usesequentialfitting = validateUseSequentialFitting(usesequentialfitting)
%usesequentialfitting = validateUseSequentialFitting(usesequentialfitting) 
%   takes a potential value for 'UseSequentialFitting' and validates it. 
%
%   What is checked?
%
%   (1) usesequentialfitting must be a logical scalar.

        usesequentialfitting = ...
            GeneralizedLinearMixedModel.validateLogicalScalar(usesequentialfitting,'stats:GeneralizedLinearMixedModel:BadUseSequentialFitting');
        
    end % end of validateUseSequentialFitting.                                  
                                          
    function ynew = validateYNew(ynew,N,subset)
%ynew = validateYNew(ynew,N,subset) takes a potential new response vector
%   ynew and returns the validated vector as a column vector. subset is a
%   logical vector of length N.
%
%   What is checked?
%
%   (1) ynew must be a numeric, real vector of length N.
%
%   (2) ynew(subset) must not have any NaN values. subset represents the
%   set of observations used in the original fit.

        % 1. Validate ynew.
        isok = isnumeric(ynew) & isreal(ynew) & isvector(ynew) & (length(ynew) == N);        
        assertThat(isok,'stats:GeneralizedLinearMixedModel:Bad_YNew',num2str(N));

        hasnans = any(isnan(ynew(subset)));
        assertThat(~hasnans,'stats:GeneralizedLinearMixedModel:Bad_YNew',num2str(N));
        
        % 2. Convert ynew to column vector.
        if size(ynew,1) == 1
            ynew = ynew';
        end
        
    end % end of validateYNew.
    
end

% Public, static fitting methods.
methods(Static, Access='public', Hidden=true)
    
    function model = fit(ds,formula,varargin)
%   The FIT method is not intended to be called directly. Use FITGLME to 
%   create a GeneralizedLinearMixedModel by fitting to data.
%
%   See also FITGLME.

        % 1. Dataset input always converted to table.
        if isa(ds,'dataset')
            ds = dataset2table(ds);
        end

        % 2. Ensure that ds and formula are sensible.
            assertThat(isa(ds,'table'),'stats:GeneralizedLinearMixedModel:Fit_firstinput');            
            assertThat(internal.stats.isString(formula),'stats:GeneralizedLinearMixedModel:Fit_secondinput');
        
        % 3. Create an empty GLME model.
        model = GeneralizedLinearMixedModel();
        
        % 4. Try to parse the formula using all variables in ds. If the
        % response vector is logical, convert it to double in ds.
        model.Formula = classreg.regr.LinearMixedFormula(formula,ds.Properties.VariableNames);        
        yresp = ds.(model.Formula.ResponseName);
        if islogical(yresp)
            ds.(model.Formula.ResponseName) = double(yresp);
        end
        
        % 5. How many grouping variable? How many total observations?
            R = length(model.Formula.RELinearFormula);
            N = size(ds,1);
            
        % 6. Initialize default values for optional parameters.
        
            % 6.1 dfltCovariancePattern
                dfltCovariancePattern = cell(R,1);
                % Default is 'Full' with Cholesky for each grouping variable.
                dfltCovariancePattern(1:R) = ...
                    {classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_FULLCHOLESKY}; 

            % 6.2 dfltFitMethod
                % Default is maximum pseudo likelihood.
                dfltFitMethod = 'mpl'; 
                
            % 6.3 dfltWeights and dfltExclude
                % Default is ones(N,1), N = size(ds,1).
                dfltWeights = ones(N,1);   
                % Default is to use all (non-NaN) obs.
                dfltExclude = false(N,1);  

            % 6.4 dfltDummyVarCoding
            dfltDummyVarCoding = 'reference';

            % 6.5 dfltOptimizer
            dfltOptimizer = 'quasinewton';
            
            % 6.6 dfltOptimizerOptions. Do not change this. We use an
            % empty value of optimizeroptions to indicate that no options
            % were supplied in the validateOptimizerAndOptions method.
            dfltOptimizerOptions = [];
            
            % 6.7 dfltStartMethod
            dfltStartMethod = 'default';
            
            % 6.8 dfltVerbose
            dfltVerbose = false;
        
            % 6.9 dfltCheckHessian
            dfltCheckHessian = false;
            
            % 6.10 dfltDistribution
            dfltDistribution = 'normal';
            
            % 6.11 dfltBinomialSize
            dfltBinomialSize = ones(N,1);
            
            % 6.12 dfltLink. We use an empty value to indicate that the
            % canonical link for the specified distribution should be used.
            dfltLink = [];
            
            % 6.13 dfltOffset
            dfltOffset = zeros(N,1);
            
            % 6.14 dfltDispersionFlag
            dfltDispersionFlag = false;
            
            % 6.15 dfltPLIterations
            dfltPLIterations = 100;
            
            % 6.16 dfltPLTolerance
            dfltPLTolerance = 1e-8;
            
            % 6.17 dfltMuStart. Empty value indicates 'MuStart' was not
            % supplied.
            dfltMuStart = [];
            
            % 6.18 dfltInitPLIterations
            dfltInitPLIterations = 10;
            
            % 6.19 dfltEBMethod
            dfltEBMethod = 'default';
            
            % 6.20 dfltEBOptions. Empty value indicates no 'EBOptions' were
            % supplied.
            dfltEBOptions = []; 
            
            % 6.21 dfltCovarianceMethod
            dfltCovarianceMethod = 'Conditional';
            
            % 6.22 dfltUseSequentialFitting
            dfltUseSequentialFitting = false;
            
        % 7. Process optional parameter name/value pairs.        
        paramNames =   {'CovariancePattern',   'FitMethod',   'Weights',   'Exclude',   'DummyVarCoding',   'Optimizer',   'OptimizerOptions',   'StartMethod',   'Verbose',   'CheckHessian',   'Distribution',   'BinomialSize',   'Link',   'Offset',   'DispersionFlag',   'PLIterations',   'PLTolerance',   'MuStart',   'InitPLIterations',   'EBMethod',   'EBOptions',   'CovarianceMethod',   'UseSequentialFitting'};
        paramDflts = {dfltCovariancePattern, dfltFitMethod, dfltWeights, dfltExclude, dfltDummyVarCoding, dfltOptimizer, dfltOptimizerOptions, dfltStartMethod, dfltVerbose, dfltCheckHessian, dfltDistribution, dfltBinomialSize, dfltLink, dfltOffset, dfltDispersionFlag, dfltPLIterations, dfltPLTolerance, dfltMuStart, dfltInitPLIterations, dfltEBMethod, dfltEBOptions, dfltCovarianceMethod, dfltUseSequentialFitting};
        [covariancepattern,fitmethod,weights,exclude,dummyvarcoding,optimizer,optimizeroptions,startmethod,verbose,checkhessian,distribution,binomialsize,linkSpec,offset,dispersionflag,pliterations,pltolerance,mustart,initpliterations,ebmethod,eboptions,covariancemethod,usesequentialfitting,setflag] = internal.stats.parseArgs(paramNames, paramDflts, varargin{:});                           
        
        % 8. Validate parameter values except covariancepattern.
        fitmethod      = model.validateFitMethod(fitmethod);        
        distribution   = model.validateDistribution(distribution);                               
        weights        = model.validateWeights(weights,N,distribution);                
        exclude        = model.validateExclude(exclude,N);
        dummyvarcoding = model.validateDummyVarCoding(dummyvarcoding);        
        [optimizer,optimizeroptions] ...
                       = model.validateOptimizerAndOptions(optimizer,optimizeroptions);                                                                             
        startmethod    = model.validateStartMethod(startmethod);                      
        verbose        = model.validateVerbose(verbose);               
        checkhessian   = model.validateCheckHessian(checkhessian);                    
        binomialsize   = model.validateBinomialSize(binomialsize,N);                              
        if isempty(linkSpec)
            % No link specified - use a default based on the specified distribution.
            linkSpec         = model.defaultLink(distribution);
        end
        linkStruct           = classreg.regr.lmeutils.StandardGeneralizedLinearMixedModel.validateLink(linkSpec,fitmethod);          
        offset               = model.validateOffset(offset,N);          
        dispersionflag       = model.validateDispersionFlag(dispersionflag);          
        pliterations         = model.validatePLIterations(pliterations);          
        pltolerance          = model.validatePLTolerance(pltolerance);         
        mustart              = model.validateMuStart(mustart,distribution,N);          
        initpliterations     = model.validateInitPLIterations(initpliterations);          
        [ebmethod,eboptions] = model.validateEBParameters(ebmethod,eboptions);          
        covariancemethod     = model.validateCovarianceMethod(covariancemethod);          
        usesequentialfitting = model.validateUseSequentialFitting(usesequentialfitting);
              
        % 9. verbose will override the Display field in optimizeroptions.
        % By default, we won't show the PL optimizer display.
        showploptimizerdisplay = false;
        if (setflag.Verbose == true)
            % User has supplied a 'Verbose' input. We also know that
            % verbose is either 0, 1 or 2.
            switch verbose
                case 0
                    optimizeroptions.Display = 'off';
                    showploptimizerdisplay = false;
                case 1
                    optimizeroptions.Display = 'iter';
                    showploptimizerdisplay = false;
                case 2
                    optimizeroptions.Display = 'iter';
                    showploptimizerdisplay = true;
            end            
        end
               
        % 10. At this point, optimizer and optimizeroptions should be
        % synchronized with each other and optimizeroptions should be
        % completely filled out with the default options overridden by the
        % user supplied options.
        
        % 11. Store parameter values in the object. weights and exclude are
        % stored in ObservationInfo and covariancepattern is stored in the
        % object after validation later.        
        model.FitMethod            = fitmethod;
        model.DummyVarCoding       = dummyvarcoding;        
        model.Optimizer            = optimizer;
        model.OptimizerOptions     = optimizeroptions;        
        model.StartMethod          = startmethod;
        model.CheckHessian         = checkhessian;        
        model.Distribution         = distribution;        
        model.BinomialSize         = binomialsize;        
        model.Link                 = linkStruct;        
        model.Offset               = offset;        
        model.DispersionFlag       = dispersionflag;        
        model.PLIterations         = pliterations;
        model.PLTolerance          = pltolerance;        
        model.MuStart              = mustart;        
        model.InitPLIterations     = initpliterations;        
        model.EBMethod             = ebmethod;        
        model.EBOptions            = eboptions;        
        model.CovarianceMethod     = covariancemethod;        
        model.UseSequentialFitting = usesequentialfitting;        
        model.ShowPLOptimizerDisplay = showploptimizerdisplay;
        
        % 12. Provisionally fill out VariableInfo and ObservationInfo. We
        % will alter PredLocs and RespLocs later in selectVariables. This
        % also indirectly validates weights and exclude arguments.
            model.PredictorTypes = 'mixed'; % Grouping vars are predictors.
            model = assignData(model,ds,[],weights,[],...
                model.Formula.VariableNames,exclude);          
        
        % 13. Select our variables and observations. After this point
        % model.ObservationInfo.Subset marks the subset of data to be used
        % for fitting after removing excluded and missing values.
            model = selectVariables(model);
            model = selectObservations(model,exclude);                       
                      
        % 14. Get y, FixedInfo, RandomInfo and GroupingInfo for the
        % observations that are going to be used in the fit. BinomialSize, 
        % Offset and MuStart are synchronized with y. Error if subset is 
        % all false.
                        subset = model.ObservationInfo.Subset;
                        if all(subset == false)
                           error(message('stats:LinearMixedModel:NoUsableObservations')); 
                        end
                      dssubset = ds(subset,:);
                       model.y = extractResponse(model,dssubset);
               model.FixedInfo = extractFixedInfo(model,dssubset);
              model.RandomInfo = extractRandomInfo(model,dssubset);
            model.GroupingInfo = extractGroupingInfo(model,dssubset);                       
            clear('dssubset');
            
            model.BinomialSize = model.BinomialSize(subset);
            model.Offset       = model.Offset(subset);
            if ~isempty(model.MuStart)
                model.MuStart = model.MuStart(subset);
            end
            
        % 15. If we have a binomial distribution, modify y by dividing it
        % by the binomialsize.
        if strcmpi(distribution,'binomial')
            model.y = model.y ./ model.BinomialSize;
        end
            
        % 16. Now that we know the size of various random effects design
        % matrices, validate covariancepattern.
        covariancepattern = GeneralizedLinearMixedModel.validateCovariancePattern...
            (covariancepattern,R,model.RandomInfo.q);

        % 17. Store covariancepattern in CovariancePattern property.
        if ( setflag.CovariancePattern == false )
            for k = 1:R
               if model.RandomInfo.q(k) == 1
                    covariancepattern{k} = classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_ISOTROPIC;
               end
            end
        end
        model.CovariancePattern = covariancepattern;                       
        
        % 18. Fit the model by calling StandardGeneralizedLinearMixedModel in fitter.        
        %   doFit is a FitObject method that does the following:
        %   (1) model = selectVariables(model);            % Set PredLocs, RespLoc and update VariableInfo.InModel
        %   (2) model = selectObservations(model,exclude); % Update ObservationInfo.Missing, .Excluded and .Subset
        %   (3) model = fitter(model);                     % Do the actual fitting.
        %   (4) model = postFit(model);                    % Do post fitting.           
        model = doFit(model);

        % 19. Omit excluded points from range.
        model = updateVarRange(model);               
        
    end % end of fit.
        
end % end of methods(Static, Access='public', Hidden).

% Public methods of LinearMixedModel.
methods(Access='public')
    
    function [D,gnames] = designMatrix(model,designtype,gnumbers)
%designMatrix Extracts the fixed or random effects design matrices.
%   X = designMatrix(GLME) or designMatrix(GLME,'Fixed') returns the N-by-P
%   fixed effects design matrix X for the generalized linear mixed effects
%   model GLME where N is the number of observations and P is the number of
%   fixed effects.
%
%   D = designMatrix(GLME,'Random') returns the overall random effects
%   design matrix corresponding to a vector B of all random effects in the
%   generalized linear mixed effects model GLME. Suppose GLME has R
%   grouping variables named g_1,...,g_R. Let Q_1,...,Q_R be the length of
%   random effects vectors associated with g_1,...,g_R respectively. Also,
%   suppose g_1,...,g_R have levels M_1,...,M_R respectively. Then B will
%   be a column vector of length Q_1*M_1 + Q_2*M_2 + ... + Q_R*M_R. B is
%   made by concatenating the empirical Bayes predictors of random effects
%   vectors corresponding to each level of each grouping variable in the
%   following order:
%
%   g_1 level 1, g_1 level 2,..., g_1 level M_1 followed by
%   g_2 level 1, g_2 level 2,..., g_2 level M_2 followed by
%      .            .                .           
%   g_R level 1, g_R level 2,..., g_R level M_R
%
%   If GLME has N observations then D will be of size N-by-length(B) such
%   that D*B is a N-by-1 vector that represents the contribution of all
%   random effects to the linear predictor of the fitted GLME model.
%
%   DSUB = designMatrix(GLME,'Random',GNUMBERS) returns a submatrix of the
%   full random effects design matrix. GNUMBERS is a length K integer array
%   with elements in the range 1 to R. DSUB is a subset of the full random
%   effects design matrix corresponding to the grouping variable names
%   indicated by integers in GNUMBERS. For example, suppose GNUMBERS is
%   [1,R] then this specifies only grouping variables g_1 and g_R. Let BSUB
%   be a vector made by concatenating empirical Bayes predictors of random 
%   effects vectors corresponding to each level of g_1 and g_R in the 
%   following order:
%
%   g_1 level 1, g_1 level 2,..., g_1 level M_1 followed by         
%   g_R level 1, g_R level 2,..., g_R level M_R
%
%   Then DSUB will be a N-by-length(BSUB) matrix such that DSUB*BSUB
%   represents the contribution of all random effects corresponding to
%   grouping variables g_1 and g_R to the linear predictor of the fitted
%   GLME model. If GNUMBERS is empty, the full random effects design matrix
%   is returned.
%
%   [DSUB,GNAMES] = designMatrix(GLME,DESIGNTYPE,GNUMBERS) also returns a 
%   K-by-1 cell array containing the names of grouping variables
%   corresponding to integers in GNUMBERS if DESIGNTYPE is 'Random'. If
%   DESIGNTYPE is 'Fixed' then GNAMES is [] and GNUMBERS is ignored.
%        
%   Example: Fit a quadratic model with continuous and categorical fixed
%            effects. Look at a few rows of the design matrix showing the
%            constant term, the Weight term, and dummy variables for the
%            Origin term.
%      load carsmall
%      T = table(MPG,Weight,Model_Year,Origin);
%      glme = fitglme(T,'MPG ~ Weight + Origin + (1|Model_Year)','Distribution','Gaussian');
%      dm = designMatrix(glme);
%      disp(dm(25:30,:))
%
%   See also RESPONSE, FITTED, RESIDUALS.

        if nargin < 2
            [D,gnames] = designMatrix@classreg.regr.LinearLikeMixedModel(model);
        elseif nargin < 3
            [D,gnames] = designMatrix@classreg.regr.LinearLikeMixedModel(model,designtype);
        else
            [D,gnames] = designMatrix@classreg.regr.LinearLikeMixedModel(model,designtype,gnumbers);
        end       
        
    end % end of designMatrix.
            
    function [beta,betanames,fetable] = fixedEffects(model,varargin)
%fixedEffects Returns estimates of fixed effects and related statistics.
%   BETA = fixedEffects(GLME) returns a vector of estimated fixed effects
%   from a fitted generalized linear mixed effects model GLME.
%
%   [BETA,BETANAMES] = fixedEffects(GLME) also returns a table array
%   BETANAMES containing the name of each fixed effects coefficient in
%   BETA.
%
%   [BETA,BETANAMES,STATS] = fixedEffects(GLME) also returns a table array
%   STATS containing the estimates of fixed effects and related statistics.
%   Table STATS has one row for each fixed effect and the following
%   columns:
%
%       Name        Name of the fixed effect coefficient
%       Estimate    Estimated coefficient value
%       SE          Standard error of the estimate
%       tStat       t statistic for a test that the coefficient is zero
%       DF          Estimated degrees of freedom for the t statistic
%       pValue      p-value for the t statistic
%       Lower       Lower limit of a 95% confidence interval
%       Upper       Upper limit of a 95% confidence interval
%
%   SE does not account for the uncertainty in estimating the covariance
%   parameters if 'CovarianceMethod' is 'Conditional' in the call to
%   FITGLME. 
%
%   [BETA,BETANAMES,STATS] = fixedEffects(GLME,'PARAM','VALUE',...) also
%   specifies optional parameter name/value pairs to control the fixed
%   effects statistics:
%  
%           'Name'     'Value'
%       'DFMethod'     Specifies the method to use for computing the
%                      approximate degrees of freedom (DF) for the t
%                      statistics that test the fixed effects coefficients
%                      against 0. Options are 'Residual' and 'None'. If
%                      'DFMethod' is 'Residual', the DF values are assumed
%                      to be constant and equal to (N-P) where N is the
%                      number of observations and P is the number of fixed
%                      effects. If 'DFMethod' is 'None', then all DF values
%                      are set to infinity. Default is 'Residual'.
%
%          'Alpha'     ALPHA, a number between 0 and 1. The columns 'Lower'
%                      and 'Upper' in table STATS will correspond to the
%                      lower and upper limits respectively of 100*(1-ALPHA)
%                      confidence intervals for fixed effects. Default is
%                      ALPHA=0.05 for 95% confidence intervals.
%
%   Example: Fit a model with a fixed effect and a random effect. Get the
%            estimated fixed effects (intercept and slope).
%      load carsmall
%      T = table(MPG,Weight,Model_Year);
%      glme = fitglme(T,'MPG ~ Weight + (1|Model_Year)')
%      fixedEffects(glme)
%
%   See also coefCI, coefTest, randomEffects.

        switch nargout
            case {0,1}
                beta                     = fixedEffects@classreg.regr.LinearLikeMixedModel(model,varargin{:});
            case 2
                [beta,betanames]         = fixedEffects@classreg.regr.LinearLikeMixedModel(model,varargin{:});
            case 3
                [beta,betanames,fetable] = fixedEffects@classreg.regr.LinearLikeMixedModel(model,varargin{:});
        end
        
    end % end of fixedEffects.

    function [b,bnames,retable] = randomEffects(model,varargin)
%randomEffects Extract estimates of random effects and related statistics.
%   B = randomEffects(GLME) returns estimates of the empirical Bayes
%   predictors (EBPs) of all random effects in the generalized linear mixed
%   effects model GLME conditional on the estimated covariance parameters
%   and the observed response. The EBPs of random effects are approximated
%   by the mode of the empirical posterior distribution of the random
%   effects given the estimated covariance parameters and the observed
%   response. Suppose GLME has R grouping variables named g_1,...,g_R. Let
%   Q_1,...,Q_R be the length of random effects vectors associated with
%   g_1,...,g_R respectively. Also, suppose g_1,...,g_R have levels
%   M_1,...,M_R respectively. Then B will be a column vector of length
%   Q_1*M_1 + Q_2*M_2 + ... + Q_R*M_R. B is made by concatenating the EBPs
%   of random effects vectors corresponding to each level of each grouping
%   variable in the following order:
%
%   g_1 level 1, g_1 level 2,..., g_1 level M_1 followed by
%   g_2 level 1, g_2 level 2,..., g_2 level M_2 followed by
%      .            .                .           
%   g_R level 1, g_R level 2,..., g_R level M_R
%
%   [B,BNAMES] = randomEffects(GLME) also returns a table array BNAMES
%   containing the names of the coefficients in B.
%
%   [B,BNAMES,STATS] = randomEffects(GLME) also returns a table array STATS
%   containing the estimates of random effects and related statistics.
%   Table STATS has one row for each random effect associated with a
%   particular grouping variable, level and predictor name with the
%   following columns:
%
%       Group       Grouping variable associated with this random effect
%       Level       Level within the grouping variable
%       Name        Name of the random effect coefficient
%       Estimate    Empirical Bayes predictor of random effect
%       SEPred      Square root of CMSEP given covariance parameters and 
%                   response
%       tStat       t statistic for a test that the random effect is zero
%       DF          Estimated degrees of freedom for the t statistic
%       pValue      p-value for the t statistic
%       Lower       Lower limit of a 95% Bayesian credible interval
%       Upper       Upper limit of a 95% Bayesian credible interval
%
%   SEPred above contains the square root of conditional mean squared error
%   of prediction (CMSEP) conditional on the estimated covariance
%   parameters and the observed response. SEPred accounts for the
%   variability in estimating the fixed effects but not the covariance
%   parameters. An alternative interpretation of SEPred is the posterior
%   standard deviation of random effects given the estimated covariance
%   parameters and the observed response.
%
%   [B,BNAMES,STATS] = randomEffects(GLME,'PARAM','VALUE',...) specifies
%   optional parameter name/value pairs to control the calculation of
%   random effects statistics:
%  
%           'Name'     'Value'
%       'DFMethod'     Specifies the method to use for computing the
%                      approximate degrees of freedom (DF) for the t
%                      statistics that test the random effects coefficients
%                      against 0. Options are 'Residual' and 'None'. If
%                      'DFMethod' is 'Residual', the DF values are assumed
%                      to be constant and equal to (N-P) where N is the
%                      number of observations and P is the number of fixed
%                      effects. If 'DFMethod' is 'None', then all DF values
%                      are set to infinity. Default is 'Residual'.
%
%          'Alpha'     ALPHA, a number between 0 and 1. The columns 'Lower'
%                      and 'Upper' in table STATS will correspond to the
%                      lower and upper limits respectively of 100*(1-ALPHA)
%                      credible intervals for random effects. Default is
%                      ALPHA=0.05 for 95% credible intervals.
%
%   Example: Fit a model with a fixed effect and random effects. Get the
%            estimated random effects. They are highly negatively
%            correlated. Verify this by plotting them.
%      load carsmall
%      T = table(MPG,Weight,Model_Year);
%      glme = fitglme(T,'MPG ~ Weight + (Weight|Model_Year)','Distribution','Normal')
%      re = randomEffects(glme)
%      plot(re(1:2:end),re(2:2:end),'rs')
%
%   See also coefCI, coefTest, fixedEffects.

%   References:  
%   (1) James G. Booth and James P. Hobert (1998), Standard Errors of
%   Prediction in Generalized Linear Mixed Models, Journal of the American
%   Statistical Association, 93(441), 262-272.

        switch nargout
            case {0,1}
                b                  = randomEffects@classreg.regr.LinearLikeMixedModel(model,varargin{:});
            case 2
                [b,bnames]         = randomEffects@classreg.regr.LinearLikeMixedModel(model,varargin{:});
            case 3
                [b,bnames,retable] = randomEffects@classreg.regr.LinearLikeMixedModel(model,varargin{:});
        end        
        
    end % end of randomEffects.

    function [PSI,mse,covtable] = covarianceParameters(model,varargin)
%covarianceParameters Extract estimated covariance parameters of the model.
%   PSI = covarianceParameters(GLME) extracts the estimated covariance
%   parameters that parameterize the prior covariance of random effects. If
%   the generalized linear mixed effects model GLME has R grouping
%   variables named g_1,...,g_R then the output PSI will be a R-by-1 cell
%   array such that PSI{i} will contain the covariance matrix of random
%   effects associated with grouping variable g_i. The order in which
%   grouping variables are assigned numbers 1 to R is the same order in
%   which grouping variables are entered into the FITGLME function.
%
%   [PSI,DISPERSION] = covarianceParameters(GLME) also extracts an estimate
%   of the dispersion parameter.
%
%   [PSI,DISPERSION,STATS] = covarianceParameters(GLME) also returns a cell
%   array STATS of length (R+1) containing covariance parameters and
%   related statistics. STATS{i} is a table array containing statistics on
%   covariance parameters for the i-th grouping variable. STATS{R+1}
%   contains statistics on the dispersion parameter. STATS{i} contains
%   columns that name each covariance parameter as well as the following
%   columns:
%
%       Group       Grouping variable name
%       Estimate    Estimate of the covariance parameter
%       Lower       Lower limit of a 95% confidence interval
%       Upper       Upper limit of a 95% confidence interval
%
%   NOTE: It is recommended that the presence or absence of covariance
%   parameters in GLME be tested using the COMPARE method, which uses a
%   likelihood ratio test. The confidence intervals in STATS are based on a
%   Laplace approximation to the log likelihood of the GLME model.
%
%   [PSI,DISPERSION,STATS] = covarianceParameters(GLME,PARAM1,VALUE1,...)
%   specifies optional name/value pairs as follows:
%
%           'Name'     'Value'
%          'Alpha'     ALPHA, a number between 0 and 1. The columns 'Lower'
%                      and 'Upper' in table STATS{i} will correspond to
%                      the lower and upper limits respectively of
%                      100*(1-ALPHA) confidence intervals for covariance
%                      parameters. Default is ALPHA = 0.05 for 95% 
%                      confidence intervals. 
%
%   Example: Fit a model with two correlated random effects, and get their
%            estimated covariance matrix.
%      load carsmall
%      T = table(MPG,Weight,Model_Year);
%      glme = fitglme(T,'MPG ~ Weight + (Weight|Model_Year)','Distribution','Normal')
%      V = covarianceParameters(glme);
%      V{1}
%
%   See also COMPARE, fixedEffects, randomEffects.

        switch nargout            
            case {0,1}
                PSI                = covarianceParameters@classreg.regr.LinearLikeMixedModel(model,varargin{:});
            case 2
                [PSI,mse]          = covarianceParameters@classreg.regr.LinearLikeMixedModel(model,varargin{:});
            case 3
                [PSI,mse,covtable] = covarianceParameters@classreg.regr.LinearLikeMixedModel(model,varargin{:});
        end

    end % end of covarianceParameters.

    function yfit = fitted(model,varargin)
%FITTED Fitted response from a generalized linear mixed effects model.
%   MUFIT = FITTED(GLME) returns the N-by-1 vector representing the fitted
%   conditional mean of the generalized linear mixed effects model GLME,
%   where N is the number of observations. If the GLME model has a N-by-P
%   fixed effects design matrix X, an estimated P-by-1 fixed effects vector
%   BETA, a N-by-Q overall random effects design matrix Z and empirical
%   Bayes predictor B of the overall random effects vector, then ETA =
%   X*BETA + Z*B + OFFSET is the linear predictor of the GLME model. If
%   ginv is the inverse of specified link function then MUFIT = ginv(ETA).
%
%   MUFIT = FITTED(GLME,PARAM1,VALUE1,...) accepts optional name/value
%   pairs as follows:
%
%           'Name'     'Value'
%     'Conditional'     Either true or false.  If false then contribution
%                       from random effects is not included. In this case,
%                       MUFIT = ginv(ETA) and ETA = X*BETA + OFFSET.
%                       Default is true.
%
%   Example: Fit a model and plot the residuals vs. fitted values, grouped by
%            Origin.
%      load carsmall
%      T = table(MPG,Weight,Model_Year);
%      glme = fitglme(T,'MPG ~ Weight + (1|Model_Year)','Distribution','Normal');
%      r = residuals(glme);
%      f = fitted(glme);
%      gscatter(f,r,Origin)
%
%   See also RESIDUALS, RESPONSE, designMatrix.

        % 1. Default parameter values.
        dfltConditional = true;

        % 2. Optional parameter names and their default values.
        paramNames = {  'Conditional'};
        paramDflts = {dfltConditional};
           
        % 3. Parse optional parameter name/value pairs.
        wantconditional = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});   

        % 4. Validate optional parameter values.
        wantconditional = model.validateConditional(wantconditional);

        % 5. Forward the call to StandardGeneralizedLinearMixedModel object slme.
        yr = fitted(model.slme,wantconditional);       
                
        % 6. Fill in NaNs in the output for observations not used in fit.
        % TO DO: LinearModel does not fill in NaNs for missing obs.
        subset = model.ObservationInfo.Subset;
        yfit = NaN(length(subset),1);
        yfit(subset) = yr;  
        
    end % end of fitted.

    function res = residuals(model,varargin)
%RESIDUALS Residuals from a fitted generalized linear mixed effects model.
%   R = RESIDUALS(GLME) returns the N-by-1 vector of raw conditional
%   residuals from a fitted generalized linear mixed effects model GLME
%   where N is the number of observations.
%
%   R = RESIDUALS(GLME,PARAM1,VALUE1,...) accepts optional name/value
%   pairs as follows:
%
%           'Name'     'Value'
%     'Conditional'     Either true or false. If false then marginal 
%                       residuals are returned. Default is true. 
%
%    'ResidualType'     Valid values are 'Raw' and 'Pearson'. Default is 
%                       'Raw'.
%
%   Suppose the GLME model has a N-by-1 response vector Y, a N-by-P fixed
%   effects design matrix X, an estimated P-by-1 fixed effects vector
%   BETA_HAT, an N-by-Q overall random effects design matrix Z and Q-by-1
%   empirical Bayes predictor B_HAT of random effects. If ginv is the
%   inverse of the specified link function, define:
%
%               RC = Y - ginv(X*BETA_HAT + Z*B_HAT + OFFSET)
%           RCTRUE = Y - ginv(  X*BETA   + Z*B     + OFFSET)
%
%               RM = Y - ginv(X*BETA_HAT + OFFSET)
%
%   Here BETA is the true fixed effects vector and B is the random effects
%   vector. The table below shows how the elements of residual vector are
%   calculated for each observation i and for various combinations of
%   values of 'Conditional' and 'ResidualType'. STD(RCTRUE(i)) is the
%   standard deviation of RCTRUE(i) given B evaluated at BETA = BETA_HAT
%   and B = B_HAT while STD(RMTRUE(i)) is the standard deviation of
%   RCTRUE(i) given B evaluated at BETA = BETA_HAT and B = 0.
%
%     ResidualType          Conditional=true         Conditional=false
%     ------------          ----------------         -----------------
%        'Raw'                   RC(i)                    RM(i)
%      'Pearson'          RC(i)/STD(RCTRUE(i))      RM(i)/STD(RMTRUE(i))
%
%   Example: Fit a model and plot the residuals vs. fitted values, grouped by
%            Origin.
%      load carsmall
%      T = table(MPG,Weight,Model_Year);
%      glme = fitglme(T,'MPG ~ Weight + (1|Model_Year)','Distribution','Normal');
%      r = residuals(glme);
%      f = fitted(glme);
%      gscatter(f,r,Origin)
%
%   See also FITTED, RESPONSE, designMatrix.

        % 1. Default parameter values.
        dfltConditional = true;
        dfltResidualType = 'Raw';

        % 2. Optional parameter names and their default values.
        paramNames = {  'Conditional',   'ResidualType'};
        paramDflts = {dfltConditional, dfltResidualType};
           
        % 3. Parse optional parameter name/value pairs.
        [wantconditional,residualtype] ...
            = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});   

        % 4. Validate optional parameter values.
        wantconditional = model.validateConditional(wantconditional);
        residualtype = model.validateResidualType(residualtype);
        
        % 5. Forward the call to member variable slme.
        r = residuals(model.slme,wantconditional,residualtype);
        
        % 6. Fill in NaNs in the output for observations not used in fit.
        subset = model.ObservationInfo.Subset;
        res = NaN(length(subset),1);
        res(subset) = r;        
        
    end % end of residuals.
    
    function table = compare(model,altmodel,varargin)
%COMPARE Compare fitted linear mixed effects models.
%   TABLE = COMPARE(GLME,ALTGLME) performs a likelihood ratio test
%   comparing generalized linear mixed effects models GLME and ALTGLME,
%   which have been fit to the same response vector but with different
%   model specifications. GLME must be nested within ALTGLME, i.e., GLME is
%   obtained from ALTGLME by setting some parameters in ALTGLME to fixed
%   values (such as 0). The output TABLE is a table array containing the
%   results of the likelihood ratio test. TABLE has 2 rows, the first row
%   is for GLME, and the second row is for ALTGLME. TABLE has the following
%   columns:
%
%           Model       The name of the model.
%           DF          The number of free parameters in the model.
%           AIC         Akaike information criterion for the model.
%           BIC         Bayesian information criterion for the model.
%           LogLik      The maximized log-likelihood for the model.
%           LRStat      Likelihood ratio test statistic for comparing 
%                       ALTGLME versus GLME.
%           deltaDF     DF for ALTGLME minus DF for GLME.
%           pValue      p-value for the likelihood ratio test.
%
%   The null and alternative hypotheses are as follows:
%
%       H0: Observed response vector was generated by model GLME. 
%       H1: Observed response vector was generated by model ALTGLME.
%
%   A p-value for the likelihood ratio test is computed by comparing the
%   observed likelihood ratio test statistic with a chi-squared reference
%   distribution with degrees of freedom deltaDF. A small p-value (e.g., 
%   < 0.05) leads to a rejection of H0 in favor of H1 and acceptance of 
%   model ALTGLME. On the other hand, a large p-value (e.g., > 0.05) 
%   reflects insufficient evidence to accept model ALTGLME.
%
%       (1) GLME and ALTGLME must be fitted using ApproximateLaplace,
%       Laplace or MPL prior to model comparison. Models fitted with REMPL
%       cannot be compared using a likelihood ratio test.
%  
%       (2) p-values computed using the chi-squared reference distribution
%       as described above can be:
%            - conservative, when testing for the presence or absence of
%              random effects terms and
%            - anti-conservative, when testing for the presence or absence
%              of fixed effects terms
%
%       Hence, testing for fixed effects should be performed either using
%       method fixedEffects or anova.
%
%   TABLE = COMPARE(GLME,ALTGLME,...,'PARAM',VALUE,...) specifies
%   additional name/value pairs for the likelihood ratio test:
%
%              'NAME'         'VALUE'
%      'CheckNesting'         Either true or false. If 'CheckNesting' is 
%                             true then we attempt to check if the smaller
%                             model GLME is nested in the bigger model
%                             ALTGLME. It is necessary that GLME be nested
%                             in ALTGLME for the theoretical likelihood
%                             ratio test to be valid. An error message is
%                             thrown if the nesting requirement is not
%                             satisfied. If 'REMPL' is used in fitting GLME
%                             and ALTGLME, then GLME and ALTGLME cannot be
%                             compared using a likelihood ratio test.
%
%   Example: Model gas mileage as a function of car weight, with a random
%            effect due to model year. Compare a model having random
%            intercepts with a fixed effects only model.
%      load carsmall
%      T = table(MPG,Weight,Model_Year);
%      glme = fitglme(T,'MPG ~ Weight','Distribution','Normal')
%      glme2 = fitglme(T,'MPG ~ Weight + (1 | Model_Year)','Distribution','Normal')
%      compare(glme,glme2)
%
%   See also fixedEffects, anova, randomEffects, covarianceParameters.

        % 1. model and altmodel are GeneralizedLinearMixedModel objects?
           model = model.validateObjectClass(   model,   'GLME','GeneralizedLinearMixedModel');
        altmodel = model.validateObjectClass(altmodel,'ALTGLME','GeneralizedLinearMixedModel');

        % 2. Default parameter values.     
        dfltCheckNesting = false;

        % 3. Optional parameter names and their default values.
        paramNames = {  'CheckNesting'};
        paramDflts = {dfltCheckNesting};
           
        % 4. Parse optional parameter name/value pairs.
        checknesting = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});   

        % 5. Validate optional parameter values.                   
        checknesting = model.validateCheckNesting(checknesting);
                    
        % 6. Get names for model and altmodel.
        modelName = inputname(1);
        altmodelName = inputname(2);
        if isempty(modelName) || isempty(altmodelName)
            modelName = 'GLME';
            altmodelName = 'ALTGLME';
        end
        
        % 7. Our expectation is that model with name modelName is the
        % smaller model and altmodel with name altmodelName is the bigger
        % model. In other words, we expect the following:
        %
        % DF_model < DF_altmodel
        %
        % model.LogLikelihood <= altmodel.LogLikelihood.
        %
        % If model and altmodel behave opposite to our expectation i.e, if
        % they satisfy:
        %
        % DF_altmodel < DF_model and
        %
        % altmodel.LogLikelihood <= model.LogLikelihood
        % 
        % then it makes sense to swap model and altmodel before proceeding
        % further. We will warn that such a swap was made.
        DF_model    = getTotalNumberOfParameters(model);
        DF_altmodel = getTotalNumberOfParameters(altmodel);
        if ( DF_altmodel < DF_model && altmodel.LogLikelihood <= model.LogLikelihood )
            
            warning(message('stats:GeneralizedLinearMixedModel:NestingCheck_modelswap',modelName,altmodelName));
            
            % Swap models.
            temp_model = model;
            model      = altmodel;
            altmodel   = temp_model;
            clear temp_model;
            
            % Swap model names.
            temp_modelName = modelName;
            modelName      = altmodelName;
            altmodelName   = temp_modelName;
            clear temp_modelName;            
        end
        
        % 8. After a possible swap, verify that altmodel attains a higher
        % maximized log likelihood compared to model.
        logliksmall = model.LogLikelihood;
        loglikbig   = altmodel.LogLikelihood;
        assertThat(loglikbig >= logliksmall,'stats:GeneralizedLinearMixedModel:NestingCheck_loglik',altmodelName,modelName);
        
        % 9. Standard LRT.                                          
            % 9.1 Ensure that model is "nested" in altmodel. This is not
            % foolproof but at least it is some protection.
            if ( checknesting == true )
                model.checkNestingRequirement(model,altmodel,modelName,altmodelName);
            end            
            % 9.2 Get the standard LRT table.
            table = model.standardLRT(model,altmodel,modelName,altmodelName);                        
        
    end % end of compare.
    
    function hout = plotResiduals(model,plottype,varargin)
% plotResiduals Plot residuals of fitted model.
%   plotResiduals(GLME,PLOTTYPE) plots the raw conditional residuals from a
%   fitted generalized linear mixed effects model GLME in a plot of type
%   PLOTTYPE. Valid values for PLOTTYPE are:
%  
%      'caseorder'     residuals vs. case (row) order
%      'fitted'        residuals vs. fitted values
%      'histogram'     histogram (default)
%      'lagged'        residual vs. lagged residual (r(t) vs. r(t-1))
%      'probability'   normal probability plot
%      'symmetry'      symmetry plot
%  
%   plotResiduals(MODEL,PLOTTYPE,'PARAM','VALUE',...) accepts optional
%   parameter name/value pairs:
%
%           'Name'     'Value'
%    'ResidualType'     Valid values are 'Raw' (default) and 'Pearson'.
%
%   For more details on the various residual types, please see help for the
%   residuals method. Both 'Raw' and 'Pearson' residuals in plotResiduals
%   are conditional residuals.
%  
%   H = plotResiduals(...) returns a handle to the lines or patches in the
%   plot.
%  
%   The PLOTTYPE argument or 'ResidualType' name-value pair can be followed
%   by parameter/value pairs to specify additional properties of the
%   primary line in the plot. For example, plotResiduals(M,'fitted','Marker','s') 
%   uses a square marker.
%     
%   For many of these plots, the data cursor tool in the figure window will
%   display the X and Y values for any data point, along with the
%   observation name or number.
%
%   Example: Fit a model and make a probability plot of the residuals.
%            There are two outliers in the upper right of the plot, which
%            you can examine using the data cursor in the plot figure.
%      load carsmall
%      obsname = strcat(num2str((1:100)'),{' '},Model);
%      T = table(MPG,Weight,Model_Year,'RowNames',obsname);
%      glme = fitglme(T,'MPG ~ Weight + (1|Model_Year)')
%      plotResiduals(glme,'probability')
%
%   See also RESIDUALS, FITTED.
        
        % 1. Prepare argument list based on number of inputs.
        if nargin < 2
            % Called like: plotResiduals(model)
            args = {}; 
        elseif nargin < 3
            % Called like: plotResiduals(model,plottype)
            args = {plottype};
        else
            % Called like: plotResiduals(model,plottype,varargin)
            args = [{plottype},varargin];
        end
        
        % 2. Call with the requested number of outputs.
        switch nargout
            case 0
                       plotResiduals@classreg.regr.LinearLikeMixedModel(model,args{:});
            case 1
                hout = plotResiduals@classreg.regr.LinearLikeMixedModel(model,args{:});
        end                
               
    end % end of plotResiduals.

    function stats = anova(model,varargin)
%ANOVA Perform hypothesis tests on fixed effect terms.
%   STATS = ANOVA(GLME) tests the significance of each fixed effect term in
%   the generalized linear mixed effects model GLME and returns a table
%   array STATS. Each fixed effect term reported in STATS is either a
%   continuous variable, a grouping variable or an interaction between two
%   or more variables (continuous or grouping). For each fixed effect term,
%   anova performs an F test (marginal test) that all coefficients
%   representing the fixed effect term are zero. If the 'CovarianceMethod'
%   is 'Conditional' in the call to FITGLME, then these F tests are
%   conditional on the estimated covariance parameters. If the
%   'CovarianceMethod' is 'JointHessian' in the call to FITGLME then these
%   F tests account for the uncertainty in estimation of covariance
%   parameters. There is one row in STATS for each fixed effect term and
%   the following columns:
%
%       Term        Name of the fixed effect term
%       FStat       F-statistic for the term
%       DF1         Numerator degrees of freedom for the F-statistic
%       DF2         Denominator degrees of freedom for the F-statistic
%       pValue      p-value for the term
%
%   To obtain tests for the Type III hypotheses, set the 'DummyVarCoding'
%   to 'effects' when fitting the model.
%
%   STATS = ANOVA(GLME,PARAM1,VALUE1,...) accepts additional parameter
%   name/value pairs:
%  
%           'Name'     'Value'
%       'DFMethod'     Specifies the method to use for computing the
%                      approximate denominator degrees of freedom DF2 for
%                      the F-statistics reported in STATS. Options are
%                      'Residual' and 'None'. If 'DFMethod' is 'Residual',
%                      the DF2 values are assumed to be constant and equal
%                      to (N-P) where N is the number of observations and P
%                      is the number of fixed effects. If 'DFMethod' is
%                      'None', then all DF2 values are set to infinity.
%                      Default is 'Residual'.
%
%   Example: Fit a model with two fixed-effect predictors and compute an
%            anova table. Note that the p-values for the intercept and
%            Weight terms are the same as in the coefficients table. The
%            p-value for the Cylinders term measures the combined
%            significance of both Cylinders coefficients.
%      load carsmall
%      T = table(MPG,Weight,Model_Year,Cylinders);
%      T.Cylinders = nominal(T.Cylinders);
%      glme = fitglme(T,'MPG ~ Weight + Cylinders + (1|Model_Year)')
%      anova(glme)
%
%   See also coefTest, coefCI, fixedEffects.

        stats = anova@classreg.regr.LinearLikeMixedModel(model,varargin{:});       

    end % end of anova.
    
    function [Y,binomsize] = response(model)
%RESPONSE Response vector for the generalized linear mixed effects model.
%   Y = RESPONSE(GLME) returns the N-by-1 response vector Y used to fit the
%   generalized linear mixed effects model GLME, where N is the number of
%   observations.
%
%   [Y,BINOMIALSIZE] = RESPONSE(GLME) also returns the binomial size
%   associated with each element of Y if the conditional distribution of
%   response given the random effects is binomial. BINOMIALSIZE is empty 
%   for other distributions.
%
%   See also FITTED, RESIDUALS, designMatrix.

        % 1. Get the subset of observations in fit.
        subset = model.ObservationInfo.Subset;
        
        % 2. Initialize Y to NaN(N,1) where N = length(subset). Then fill
        % in the response values at the observations used in the fit.
        N = length(subset);
        Y = NaN(N,1);
        Y(subset) = model.y;                

        % 3. Return binomial size if required.
        if nargout > 1
            if strcmpi(model.Distribution,'binomial')
                binomsize = NaN(N,1);
                binomsize(subset) = model.BinomialSize; 
            else
                binomsize = [];
            end
        end
        
    end % end of response.
               
    function model = refit(model,ynew)
%REFIT Refit a generalized linear mixed effects model to new response.
%   GLME = REFIT(GLME,YNEW) takes a fitted generalized linear mixed effects
%   model GLME, a new response vector YNEW of size N-by-1 (where N is the
%   number of observations used to fit GLME previously), refits the same
%   model to YNEW, and returns the updated model.
%
%   Example: Fit a model. Simulate a new response vector from the fitted 
%            model and refit the model to this new response.
%      load carsmall
%      T = table(MPG,Weight,Model_Year);
%      glme = fitglme(T,'MPG ~ Weight + (1|Model_Year)');
%      rng(0,'twister');
%      ynew = random(glme);
%      glme = refit(glme,ynew);
%      disp(glme);
%
%   See also FITTED, RESIDUALS, designMatrix.
           
        % 1. How many observations in original dataset?
        subset = model.ObservationInfo.Subset;
        N = length(subset);
        
        % 2. What distribution do we have?
        distribution = model.Distribution;
        
        % 3. Validate ynew.
        ynew = model.validateYNew(ynew,N,subset);                
        
        % 4. Reset the value of response variable in the model. For
        % 'binomial' distribution, the "fractions" in ynew need to be
        % multiplied by the original binomial size.
        if strcmpi(distribution,'binomial')            
            ynew = ynew .* model.ObservationInfo.BinomSize;
        end
        ds = model.Variables;
        ds.(model.ResponseName) = ynew;
                
        % 5. Provisionally fill out VariableInfo and ObservationInfo. We
        % will alter PredLocs and RespLocs later in selectVariables. This
        % also indirectly validates weights and exclude arguments.
        weights = model.ObservationInfo.Weights;
        exclude = model.ObservationInfo.Excluded;
            model.PredictorTypes = 'mixed'; % Grouping vars are predictors.
            model = assignData(model,ds,[],weights,[],...
                model.Formula.VariableNames,exclude);          
        
        % 6. Select our variables and observations. After this point
        % model.ObservationInfo.Subset marks the subset of data to be used
        % for fitting after removing excluded and missing values.
            model = selectVariables(model);
            model = selectObservations(model,exclude);                       
                      
        % 7. Get y for the observations that are used in the fit.                                
                    if all(subset == false)
                        error(message('stats:GeneralizedLinearMixedModel:NoUsableObservations'));
                    end
                    dssubset = ds(subset,:);
                    model.y = extractResponse(model,dssubset);                        
                    clear('dssubset');
                        
        % 8. If we have a binomial distribution, modify y by dividing it
        % by the binomialsize.
        if strcmpi(distribution,'binomial')
            model.y = model.y ./ model.BinomialSize;
        end                   
        
        % 9. Fit the model by calling StandardGeneralizedLinearMixedModel in fitter.        
        %   doFit is a FitObject method that does the following:
        %   (1) model = selectVariables(model);            % Set PredLocs, RespLoc and update VariableInfo.InModel
        %   (2) model = selectObservations(model,exclude); % Update ObservationInfo.Missing, .Excluded and .Subset
        %   (3) model = fitter(model);                     % Do the actual fitting.
        %   (4) model = postFit(model);                    % Do post fitting.           
        model = doFit(model);

        % 10. Omit excluded points from range.
        model = updateVarRange(model);         
        
    end % end of refit.
        
end % end of methods(Access='public').

end

%=== helper method for verifying assertions
function assertThat(condition,msgID,varargin)
%   assertThat(condition,msgID,varargin) takes a variable condition that is
%   either true or false, a message catalog ID msgID and optional arguments
%   required to create a message object from msgID. If condition is false,
%   an error message represented by msgID is thrown.

    if ~condition        
        % (1) Create a message object from msgID and varargin.
        try
            msg = message(msgID,varargin{:});
        catch
            % <entry key="BadMsgID">Invalid message ID: {0}.</entry>      
            error(message('stats:LinearMixedModel:BadMsgID',msgID));
        end        
        % (2) Create and throw an MException.
        ME = MException( msg.Identifier, getString(msg) );
        throwAsCaller(ME);        
    end

end % end of assertThat.

