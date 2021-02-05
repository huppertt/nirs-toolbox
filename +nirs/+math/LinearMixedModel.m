classdef (Sealed = true) LinearMixedModel < classreg.regr.LinearLikeMixedModel
%LinearMixedModel Fitted linear mixed effects model.
%   LME = FITLME(...) fits a linear mixed effects model to data. The fitted
%   model LME is a LinearMixedModel modeling a response variable as a
%   linear function of fixed effect predictors and random effect
%   predictors.
%
%   LinearMixedModel methods:
%       coefCI - Coefficient confidence intervals
%       coefTest - Linear hypothesis test on coefficients
%       predict - Compute predicted values given predictor values
%       random - Generate random response values given predictor values
%       plotResiduals - Plot of residuals
%       designMatrix - Fixed and random effects design matrices
%       fixedEffects - Stats on fixed effects
%       randomEffects - Stats on random effects
%       covarianceParameters - Stats on covariance parameters
%       fitted - Fitted response
%       residuals - Various types of residuals
%       response - Response used to fit the model
%       compare - Compare fitted models
%       anova - Marginal tests for fixed effect terms           
%       disp - Display a fitted model                                                                                                                 
%
%   LinearMixedModel properties:
%       FitMethod - Method used for fitting (either 'ML' or 'REML')
%       MSE - Mean squared error (estimate of residual variance)
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
%   See also FITLME, FITLMEMATRIX, LinearModel, GeneralizedLinearModel, NonLinearModel.    
     
%   Copyright 2012-2017 The MathWorks, Inc.

% Properties from Predictor.
properties(Dependent,GetAccess='public',Hidden=true,SetAccess='protected')    
%Fitted - Fitted (predicted) values.
%   The Fitted property is a vector of conditional fitted values. They
%   include contributions from both fixed and random effects.
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
%        'Raw'          The Raw conditional residual
%        'Pearson'      The Pearson conditional residual
%        'Standardized' The Standardized conditional residual
%
%   The "conditional" residuals are based on fitted values that include
%   contributions from both fixed and random effects.
%
%   For more details on various residuals, see the residuals method.
%
%   See also RESIDUALS.
    Residuals    

end
  
% New public properties defined by LinearMixedModel.
properties(GetAccess='public')
%FitMethod - Method used to fit the linear mixed effects model.
%   The FitMethod property is the string 'ML' if the model was fit using
%   maximum likelihood (ML) or 'REML' if the model was fit using
%   restricted maximum likelihood (REML).
%
%   See also FITLME, FITLMEMATRIX.
    FitMethod
    pmax;
%MSE - Estimated residual variance.
%   The MSE property is an estimate of the residual variance or variance of
%   the observation error term of the linear mixed effects model. MSE
%   contains either a maximum likelihood (ML) or restricted maximum
%   likelihood (REML) based estimate of residual variance depending on
%   whether ML or REML was used in fitting the model.
%
%   See also covarianceParameters.  
    MSE  
end

% New public hidden properties defined by LinearMixedModel.
properties(GetAccess='public', Hidden=true)    
%Response - Vector of response values used to fit the model
    Response    
end

% Internal constants.
properties(Constant=true,Hidden=true)
                
    AllowedFitMethods = {'ML','REML'};
        
    AllowedDFMethods = {'None','Residual','Satterthwaite'};
    
    AllowedResidualTypes = {'Raw','Pearson','Standardized'};

    AllowedOptimizers = {'fminunc','quasinewton'};
        
    AllowedStartMethods = {'random','default'};
    
end

% Private properties.
properties(Access={?classreg.regr.LinearLikeMixedModel})
%   The following variables only include observations used in the fit 
%   corresponding to ObservationInfo.Subset.        

%XYZGNames - Structure containing the original names of predictors when using the matrix version of fit.
%   XYZGNames has the following fields:
%
%   XNames = P-by-1 cell array of strings containing names of columns of X.
%            P is the number of columns in X.
%
%    YName = String containing name for Y.
%
%   ZNames = R-by-1 cell array. ZNames{k} is a q(k)-by-1 cell array of 
%            strings containing names for columns of Z{k}.
%
%   GNames = R-by-1 cell array of strings. GNames{k} contains the name of
%            grouping variable G{k}.    
    XYZGNames

%CovariancePattern - A cell array of covariance patterns.
%   If the LME model has R grouping variables, CovariancePattern is a cell
%   array of length R containing the covariance pattern used for each 
%   grouping variable. 
    CovariancePattern
    
%DummyVarCoding - The dummy variable coding used for creating design matrices.
    DummyVarCoding
    
%Optimizer - A string containing the optimizer to use for optimization.
    Optimizer
    
%OptimizerOptions - A structure or object specifying the optimization options.
    OptimizerOptions
    
%StartMethod - A string specifying how to initialize parameters for optimization.
    StartMethod    
    
%CheckHessian - A logical scalar indicating whether Hessian checks should be performed after fitting the model.    
    CheckHessian
                      
%CovarianceTable - A cell array of tables returned by the method
%   covarianceParameters. We store this info here so that we don't have to
%   recompute this in disp.
    CovarianceTable
end

% Abstract public methods from FitObject.
methods(Access='public', Hidden=true)
    
    function t = title(model)
        strLHS = model.ResponseName;
        strFunArgs = internal.stats.strCollapse(model.Formula.PredictorNames,',');
        t = sprintf( '%s = lme(%s)',strLHS,strFunArgs);
    end % end of title.        
    
    function val = feval(model,varargin) %#ok<INUSD>
        warning(message('stats:LinearMixedModel:NoFevalMethod'));
        val = [];
    end % end of feval.
        
end

% Abstract public methods from FitObject.
methods(Access='public')
        
    function disp(model)
%DISP Display a LinearMixedModel.
%   DISP(LME) displays the LinearMixedModel LME.
%
%   See also LinearMixedModel, FITLME, FITLMEMATRIX.

        % Can't display an empty model.
        if isempty(model.ObservationInfo)
            displayFormula(model);            
            error(message('stats:LinearMixedModel:NoConstructor'));
        end
        
        % (1) Display headline.
        displayHeadLine(model);
        
        % (2) Model information.
        displayModelInfo(model);
        
        % (3) Display the formula.
        displayFormula(model);

        % (4) Model fit stats.
        displayModelFitStats(model);        
            
        % (5) Fixed effect stats.
        displayFixedStats(model)            
            
        % (6) Random effects covariance parameter stats.
        displayCovarianceStats(model);               
        
    end % end of disp.
        
end

% Abstract public methods from Predictor.
methods(Access='public')

    function [ypred,yci,df] = predict(model,varargin)
%PREDICT Compute predicted values given predictor values.
%   YPRED = PREDICT(LME) computes a vector YPRED of conditional predictions
%   from the LinearMixedModel LME at the original predictors used to create
%   LME. Conditional predictions include contributions from both fixed and
%   random effects.
%
%   If LME was created using a dataset/table input via FITLME, supply a new
%   dataset/table DS to predict like this:
%
%   YPRED = PREDICT(LME,DS) uses predictor variables from the dataset/table DS.
%   DS must contain all of the predictor variables used to create LME.
%
%   If LME was created using a matrix input via FITLMEMATRIX, supply a new
%   matrix X, new matrix or cell array of matrices Z and new grouping 
%   variable or cell array of grouping variables G in one of the following
%   two forms:
%
%   YPRED = PREDICT(LME,X,Z) uses the fixed and random effects design
%   matrices X and Z respectively. Z may also be a cell array of random
%   effects design matrices.
%
%   YPRED = PREDICT(LME,X,Z,G) uses the fixed and random effects design
%   matrices X and Z respectively and the grouping variable G. Both Z and G
%   can be cell arrays of matrices and grouping variables respectively.
%
%   [YPRED,YCI] = PREDICT(...) also returns the two-column matrix YCI
%   containing 95% pointwise confidence intervals for the predicted values.
%   The lower limits of the bounds are in column 1, and the upper limits 
%   are in column 2.
%
%   [YPRED,YCI,DF] = PREDICT(...) also returns a vector DF containing the
%   degrees of freedom values used in computing the confidence intervals if
%   'Simultaneous' is false. If 'Simultaneous' is true then DF is a scalar
%   containing the DF value used for computing the simultaneous confidence
%   intervals.
%  
%   [...] = PREDICT(...,PARAM1,VALUE1,...) specifies one or more of the
%   following name/value pairs:
%
%     'Conditional'     Either true or false. If false then returns 
%                       marginal predictions which include contribution
%                       only from the fixed effects. When making
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
%      'Prediction'     Either 'curve' or 'observation'. When 'Prediction' 
%                       is 'observation', the variability due to 
%                       observation error for the new observations is
%                       included in the confidence interval calculation.
%                       Default is 'curve'.
%
%       'DFMethod'      Specifies the method to use for computing the
%                       approximate denominator degrees of freedom (DF) 
%                       when computing the confidence intervals. Options
%                       are 'Satterthwaite', 'Residual' and 'None'. If
%                       'DFMethod' is 'Satterthwaite', a Satterthwaite
%                       approximation is used to compute DF. If 'DFMethod'
%                       is 'Residual', the DF values are assumed to be
%                       constant and equal to (N-P) where N is the number
%                       of observations and P is the number of fixed
%                       effects. If 'DFMethod' is 'None', then all DF
%                       values are set to infinity. Default is 'Residual'.
%  
%          'Alpha'      A value between 0 and 1 to specify the confidence
%                       level as 100(1-ALPHA)%.  Default is 0.05 for 95%
%                       confidence.
%
%   Example: Model gas mileage as a function of car weight, with a random
%            effect due to model year.
%      load carsmall
%      ds = dataset(MPG,Weight,Model_Year);
%      lme = fitlme(ds,'MPG ~ Weight + (1|Model_Year)')
%
%      % Plot predicted values conditional on each year.
%      gscatter(ds.Weight,ds.MPG,ds.Model_Year)
%      ds2 = dataset;
%      ds2.Weight = linspace(1500,5000)';
%      ds2.Model_Year = repmat(70,100,1);
%      line(ds2.Weight,predict(lme,ds2),'color','r');
%      ds2.Model_Year(:) = 76;
%      line(ds2.Weight,predict(lme,ds2),'color','g');
%      ds2.Model_Year(:) = 82;
%      line(ds2.Weight,predict(lme,ds2),'color','b');
%
%   See also FITTED, RANDOM.
        [varargin{:}] = convertStringsToChars(varargin{:});
        if nargin < 2 || internal.stats.isString(varargin{1})            
            % (1) If no new design points are given, get conditional
            % predictions from the fitted model at the original design
            % points.
            haveDataset = true;
            ds = model.Variables;
            X = [];
            Z = [];
            G = [];
            otherArgs = varargin;     
        else                        
            % (2) Figure out whether we have a dataset/table input or matrix form
            % input.
            [haveDataset,ds,X,Z,G,otherArgs] ...
                = LinearMixedModel.handleDatasetOrMatrixInput(varargin{:});
        end
        
        % (3) Parse and validate optional input args.
            % (1) Default parameter values.
            dfltConditional = true;
            dfltSimultaneous = false;
            dfltPrediction = 'curve';            
            dfltDFMethod = 'Residual';
               dfltAlpha = 0.05;

            % (2) Optional parameter names and their default values.
            paramNames =   {'Conditional',   'Simultaneous',   'Prediction',   'DFMethod',   'Alpha'};
            paramDflts = {dfltConditional, dfltSimultaneous, dfltPrediction, dfltDFMethod, dfltAlpha};

            % (3) Parse optional parameter name/value pairs.
            [conditional,simultaneous,prediction,dfmethod,alpha] ...
                = internal.stats.parseArgs(paramNames,paramDflts,otherArgs{:});   

            % (4) Validate optional parameter values.
            conditional  = LinearMixedModel.validateConditional(conditional);
            simultaneous = LinearMixedModel.validateSimultaneous(simultaneous);
            prediction   = LinearMixedModel.validatePrediction(prediction);
            dfmethod     = LinearMixedModel.validateDFMethod(dfmethod);
            alpha        = LinearMixedModel.validateAlpha(alpha);
  
        % (4) Extract size info from model and the number of rows in X.
              p = size(model.FixedInfo.X,2);
              q = model.RandomInfo.q;
              R = model.GroupingInfo.R;
            lev = model.GroupingInfo.lev;
            if haveDataset
                M = size(ds,1);
            else
                M = size(X,1);            
            end
            
        % (5) If we have a dataset/table, validate ds. If we don't have a 
        % dataset/table, validate X, Z and G.  
        if haveDataset == true
            % (1) Validate input ds.
            predNames = model.PredictorNames;
            varNames = model.Variables.Properties.VariableNames;
            [~,predLocs] = ismember(predNames,varNames);
            dsref = model.Variables(:,predLocs);
            ds = LinearMixedModel.validateDataset(ds,'DS',dsref);                       
        else
            % (1) Ensure X, Z and G are consistent with each other and set 
            % sensible values for empty Z and/or G.
                [X,Z,G] = LinearMixedModel.validateXZG(X,Z,G,'X','Z','G');
            % (2) Validate X against original X.            
                X = LinearMixedModel.validateMatrix(X,'X',[],p);
            % (3) Validate Z against original Z.                
                Z = LinearMixedModel.validateCellVector(Z,'Z',R);
                for k = 1:R
                    ZNamek = ['Z{',num2str(k),'}'];
                    Z{k} = LinearMixedModel.validateMatrix(Z{k},ZNamek,M,q(k));
                end
            % (4) Validate G against original G.
                G = LinearMixedModel.validateCellVector(G,'G',R);
                for k = 1:R
                    GNamek = ['G{',num2str(k),'}'];
                    G{k} = LinearMixedModel.validateGroupingVar(G{k},GNamek,M);                    
                end
            % (5) Convert X, Y, Z, G into a table and use this to build
            % X and Z. Why are we doing this? Don't we already have X and
            % Z? Yes, we do but the column order in user suplied Z may be
            % different from the one used internally. So we follow the same
            % process as in the .fitmatrix method i.e., convert X, Y, Z and
            % G to a table first and then use extractFixedInfo and
            % extractRandomInfo methods.
                % (1) Original predictor names for .fitmatrix
                fepredictors = model.XYZGNames.XNames;
                     respvar = model.XYZGNames.YName;
                repredictors = model.XYZGNames.ZNames;
                    regroups = model.XYZGNames.GNames;
                % (2) Create a dummy Y for the purposes of getting the
                % table.
                Y = NaN(M,1);
                % (3) Convert X, Y, Z, G to table ds (same technique as
                % used in .fitmatrix). Then remove Y from ds.
                [ds,~] = LinearMixedModel.convertXYZGToDataset(X,Y,Z,G,fepredictors,respvar,repredictors,regroups);
                ds.(respvar) = [];                
        end
        
        % (6) Get X and Z.
            % (1) Get X.
            finfo = extractFixedInfo(model,ds);    
                X = finfo.X;
            % (2) Get Z.
            rinfo = extractRandomInfo(model,ds);
                Z = rinfo.Z; 
                
        % (7) Get cell arrays Gid and GidLevelNames for the supplied
        % observations.
        ginfo = extractGroupingInfo(model,ds);
        Gid = ginfo.Gid;
        GidLevelNames = ginfo.GidLevelNames;       
        
        % (8) Modify Gid to newGid such that in newGid{i}, ID j is mapped
        % to model.GroupingInfo.GidLevelNames{i}{j}. In other words, if we
        % associated ID 5 with level name 'Green', we must associate ID 5
        % with level 'Green' in the supplied observations as well. New
        % levels in supplied data, not previously seen will be assigned a
        % group ID of NaN.
        newGid = cell(R,1);
        for k = 1:R           
            newGid{k} = LinearMixedModel.reorderGroupIDs(Gid{k},...
                GidLevelNames{k},model.GroupingInfo.GidLevelNames{k});
        end
        
        % (9) X already contains the fixed effects design matrix. Get the
        % full sparse random effects design matrix Zs. This matrix will
        % have sum_{i=1 to R} (q(i)*lev(i)) columns and M rows.
        Zs = LinearMixedModel.makeSparseZ(Z,q,lev,newGid,M);
  
        % (10) Get ypred and yci.        
            % (1) Prediction options.                     
            %   alpha = 0.05;
            %   dfmethod = 'residual';
                wantConditional = conditional;
                if simultaneous == true
                    wantPointwise = false;
                else
                    wantPointwise = true;
                end
                if strcmpi(prediction,'curve')
                    wantCurve = true; 
                else
                    wantCurve = false; 
                end
            % (2) Call predict on model.slme.
            hasIntercept = model.Formula.FELinearFormula.HasIntercept;
            args = {X,Zs,alpha,dfmethod,...
                wantConditional,wantPointwise,wantCurve,hasIntercept};
            switch nargout
                case {0,1}
                    ypred          = predict(model.slme,args{:});
                case 2
                    [ypred,yci]    = predict(model.slme,args{:});
                case 3
                    [ypred,yci,df] = predict(model.slme,args{:});
                    % (3) Warn about 0 degrees of freedom from 'Satterthwaite'.
                    if any(df == 0) && strcmpi(dfmethod,'Satterthwaite')
                        warning(message('stats:LinearMixedModel:BadSatterthwaiteDF'));
                    end
            end
            %[ypred,yci,df] = predict(model.slme,X,Zs,alpha,dfmethod,...
            %    wantConditional,wantPointwise,wantCurve,hasIntercept);
            
    end % end of predict.
    
    function ynew = random(model,varargin)
%RANDOM Generate random response values given predictor values.
%   YNEW = RANDOM(LME) generates a vector YNEW of random values from the
%   fitted linear mixed effects model LME at the original design points.
%   YNEW is created by starting from the fixed effects part of LME and
%   adding new realizations of the random effects and observation error
%   using the fitted LME. The effect of supplied observation weights when 
%   creating the fit (if any) is taken into account.
%
%   If LME was created using a dataset/table input via FITLME, supply a new
%   dataset/table DS to random like this:
%
%   YNEW = RANDOM(LME,DS) uses predictor variables from the dataset/table DS. DS
%   must contain all of the predictor variables used to create LME.
%
%   If LME was created using a matrix input via FITLMEMATRIX, supply new
%   matrix X, new matrix or cell array of matrices Z and new grouping 
%   variable or cell array of grouping variables G in one of the following
%   two forms:
%
%   YNEW = RANDOM(LME,X,Z) uses the fixed and random effects design
%   matrices X and Z respectively. Z may also be a cell array of random
%   effects design matrices.
%
%   YNEW = RANDOM(LME,X,Z,G) uses the fixed and random effects design
%   matrices X and Z respectively and the grouping variable G. Both Z and G
%   can be cell arrays of matrices and grouping variables respectively.
%
%   Example: Fit a model. Simulate new random values for the first
%            observation.
%      load carsmall
%      ds = dataset(MPG,Weight,Model_Year);
%      lme = fitlme(ds,'MPG ~ Weight + (1|Model_Year)');
%      newds = ds(ones(10000,1),:);
%      r = random(lme,newds);
%      hist(r)
%
%      % These random values share a common new random effect due to
%      % Model_Year, but their variance is comparable to the model value.
%      v1 = var(r)
%      v2 = lme.MSE
%
%   See also PREDICT, FITTED.

        % (1) If no new design points are given, generate data from the
        % fitted model at the original design points. Include the effect of
        % originally supplied observation weights.
        if nargin < 2            
            % (1) First get W^(1/2)*ysim and then get ysim from that.
            ysim = random(model.slme,[],model.slme.X,model.slme.Z);
            w = getCombinedWeights(model,true);
            ysim = ysim ./ sqrt(w);
            
            % (2) Put NaN's into ynew for observations not used in the fit.
            subset = model.ObservationInfo.Subset;
            ynew = NaN(length(subset),1);       
            ynew(subset) = ysim;            
            
            return;            
        end

        % (2) Figure out whether we have a dataset/table input or matrix form
        % input. We may use otherArgs in the future.        
        [haveDataset,ds,X,Z,G,otherArgs] ...
            = LinearMixedModel.handleDatasetOrMatrixInput(varargin{:}); %#ok<NASGU>

        % (3) Extract size info from model and the number of rows in X.
              p = size(model.FixedInfo.X,2);
              q = model.RandomInfo.q;
              R = model.GroupingInfo.R;
              if haveDataset
                  M = size(ds,1);
              else
                  M = size(X,1);
              end
                            
        % (4) If we have a dataset/table, validate ds. If we don't have a 
        % dataset/table, validate X, Z and G. 
        if haveDataset == true
            % (1) Validate input ds.
            predNames = model.PredictorNames;
            varNames = model.Variables.Properties.VariableNames;
            [~,predLocs] = ismember(predNames,varNames);
            dsref = model.Variables(:,predLocs);
            ds = LinearMixedModel.validateDataset(ds,'DS',dsref);             
        else
            % (1) Ensure X, Z and G are consistent with each other and set 
            % sensible values for empty Z and/or G.
                [X,Z,G] = LinearMixedModel.validateXZG(X,Z,G,'X','Z','G');
            % (2) Validate X against original X.                            
                X = LinearMixedModel.validateMatrix(X,'X',[],p);
            % (3) Validate Z against original Z.                
                Z = LinearMixedModel.validateCellVector(Z,'Z',R);
                for k = 1:R
                    ZNamek = ['Z{',num2str(k),'}'];
                    Z{k} = LinearMixedModel.validateMatrix(Z{k},ZNamek,M,q(k));
                end
            % (4) Validate G against original G.
                G = LinearMixedModel.validateCellVector(G,'G',R);
                for k = 1:R
                    GNamek = ['G{',num2str(k),'}'];
                    G{k} = LinearMixedModel.validateGroupingVar(G{k},GNamek,M);                    
                end
            % (5) Convert X, Y, Z, G into a table and use this to build
            % X and Z. Why are we doing this? Don't we already have X and
            % Z? Yes, we do but the column order in user suplied Z may be
            % different from the one used internally. So we follow the same
            % process as in the .fitmatrix method i.e., convert X, Y, Z and
            % G to a table first and then use extractFixedInfo and
            % extractRandomInfo methods.
                % (1) Original predictor names for .fitmatrix
                fepredictors = model.XYZGNames.XNames;
                     respvar = model.XYZGNames.YName;
                repredictors = model.XYZGNames.ZNames;
                    regroups = model.XYZGNames.GNames;
                % (2) Create a dummy Y for the purposes of getting the
                % table.
                Y = NaN(M,1);
                % (3) Convert X, Y, Z, G to table ds (same technique as
                % used in .fitmatrix). Then delete the response variable
                % from ds.
                [ds,~] = LinearMixedModel.convertXYZGToDataset(X,Y,Z,G,fepredictors,respvar,repredictors,regroups);
                ds.(respvar) = [];                
        end                
        
        % (5) Get X and Z from ds.
            % (1) Get X.
            finfo = extractFixedInfo(model,ds);    
                X = finfo.X;
            % (2) Get Z.
            rinfo = extractRandomInfo(model,ds);
                Z = rinfo.Z;  
                
        % (6) Get a R-by-1 cell vector Gid such that Gid{i} is the group
        % IDs for grouping variable i. Also get a R-by-1 vector lev such
        % that lev(i) is the number of levels of grouping variable i in the
        % input data. Extract grouping variable info from ds. In ginfo,
        % GidLevelNames may contain levels not previously seen during model
        % fitting. That's okay, every level gets its own random effects
        % vector.
        ginfo = extractGroupingInfo(model,ds);
        Gid = ginfo.Gid;
        lev = ginfo.lev;        
            
        % (7) X already contains the fixed effects design matrix. Get the
        % full sparse random effects design matrix Zs. This matrix will
        % have sum_{i=1 to R} (q(i)*lev(i)) columns and M rows.
        Zs = LinearMixedModel.makeSparseZ(Z,q,lev,Gid,M);        
        
        % (8) Get the sum_{i=1 to R} (q(i)*lev(i))-by-1 vector bsim.
        if isempty(lev)
            bsim = zeros(0,1);
        else
            bsim = randomb(model.slme,[],lev);
        end
        
        % (9) Get epsilonsim, the residual error realization.
        epsilonsim = model.slme.sigmaHat*randn(M,1);
        
        % (10) Form the output.
        if isempty(bsim)
            ynew = X * model.slme.betaHat + epsilonsim;
        else
            ynew = X * model.slme.betaHat + Zs * bsim + epsilonsim;
        end        

    end % end of random.
    
end

% Abstract public methods from ParametricRegression.
methods(Access='public', Hidden=true)
    
    function v = varianceParam(model) % a consistent name for the dispersion/MSE/whatever
        v = model.MSE;
    end
    
end

% Public methods from ParametricRegression.
methods(Access='public')
   
    function [feci,reci] = coefCI(model,varargin)
%coefCI Confidence intervals for coefficients.
%   FECI = coefCI(LME) computes 95% confidence intervals for the fixed
%   effects parameters in the linear mixed effects model LME. The output
%   FECI is a P-by-2 matrix where P is the number of fixed effects
%   parameters in LME. Rows of FECI from top to bottom correspond
%   respectively to the P-by-1 fixed effects vector BETA displayed from top
%   to bottom in the tabular display from the fixedEffects method. Column 1
%   of FECI displays lower confidence limits and column 2 of FECI displays
%   upper confidence limits.
%
%   [FECI,RECI] = coefCI(LME) also returns 95% confidence intervals for
%   random effects parameters in LME. The output RECI is a Q-by-2 matrix
%   where Q is the total number of random effects parameters in LME. Rows
%   of RECI from top to bottom correspond respectively to the Q-by-1 random
%   effects vector B displayed from top to bottom in the tabular display
%   from the randomEffects method. Column 1 of RECI displays lower
%   confidence limits and column 2 of RECI displays upper confidence
%   limits.
%
%   [FECI,RECI] = coefCI(LME,'PARAM','VALUE',...) accepts optional
%   name/value pairs:
%      
%           'Name'     'Value'
%          'Alpha'     ALPHA, a number between 0 and 1. Computes 
%                      100*(1-ALPHA)% confidence intervals. Default is
%                      ALPHA=0.05 for 95% confidence intervals.
%  
%       'DFMethod'     Specifies the method to use for computing the
%                      approximate degrees of freedom DF for confidence
%                      interval calculation. Options are 'Satterthwaite',
%                      'Residual' and 'None'. If 'DFMethod' is
%                      'Satterthwaite', a Satterthwaite approximation is
%                      used to compute DF. If 'DFMethod' is 'Residual', the
%                      DF values are assumed to be constant and equal to
%                      (N-P) where N is the number of observations and P is
%                      the number of fixed effects. If 'DFMethod' is
%                      'None', then all DF values are set to infinity.
%                      Default is 'Residual'.
%
%   Example: Fit a model with Weight as a fixed effect, and a random effect
%            due to Model_Year. Compute confidence intervals for the
%            intercept and the coefficient of Weight.
%      load carsmall
%      ds = dataset(MPG,Weight,Model_Year);
%      lme = fitlme(ds,'MPG ~ Weight + (1|Model_Year)')
%      coefCI(lme)
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
%   PVAL = coefTest(LME) computes the p-value for an F test that all fixed
%   effects parameters in the linear mixed effects model LME except the
%   intercept are zero.
%
%   PVAL = coefTest(LME,H) computes the p-value for an F test on the fixed
%   effects part of LME using the M-by-P matrix H where P is the number of
%   fixed effects parameters in LME. Each row of H represents 1 contrast
%   and the columns of H from left to right correspond respectively to the
%   P-by-1 fixed effects vector BETA displayed from top to bottom in the
%   tabular display from the fixedEffects method. The output PVAL is the
%   p-value for an F test that H*BETA = 0. To include contrasts that
%   involve the random effects, use the 'REContrast' parameter.
%
%   PVAL = coefTest(LME,H,C) also specifies a M-by-1 vector C for testing 
%   the hypothesis H*BETA = C.
%
%   PVAL = coefTest(LME,H,C,PARAM1,VALUE1,...) accepts optional name/value
%   pairs to control the calculation of PVAL.
%
%             Name     Value
%       'DFMethod'     A string that specifies the method to use for
%                      computing the approximate denominator degrees of
%                      freedom DF for the F test. If 'DFMethod' is
%                      'Satterthwaite', a Satterthwaite approximation is
%                      used to compute DF. If 'DFMethod' is 'Residual', the
%                      DF value is assumed to be equal to (N-P) where N is
%                      the number of observations and P is the number of
%                      fixed effects. If 'DFMethod' is 'None', then the DF
%                      value is taken to be infinity. Default is
%                      'Residual'.
%
%       'REContrast'   A M-by-Q matrix K where Q is the number of random 
%                      effects parameters in LME. Each row of K represents
%                      1 contrast and the columns of K from left to right
%                      correspond respectively to the Q-by-1 random effects
%                      vector B displayed from top to bottom in the tabular
%                      display from the randomEffects method. The output
%                      PVAL is the p-value for an F test H*BETA + K*B = C.
%
%   [PVAL,F,DF1,DF2] = coefTest(...) also returns the F-statistic F, the
%   numerator degrees of freedom DF1 for F, and the denominator degrees of
%   freedom DF2 for F. DF1 is equal to the number of linearly independent
%   rows in H, or [H,K] depending on the call to coefTest. The value of DF2
%   depends on the option selected for 'DFMethod'.
%
%   Example: Fit a model with two fixed-effect predictors and a random
%            effect. Test for the significance of the Cylinders term. The
%            p-value is the same as shown in the anova table.
%      load carsmall
%      ds = dataset(MPG,Weight,Model_Year,Cylinders);
%      ds.Cylinders = nominal(ds.Cylinders);
%      lme = fitlme(ds,'MPG ~ Weight + Cylinders + (1|Model_Year)')
%      p = coefTest(lme,[0 0 1 0;0 0 0 1])
%      anova(lme)
%
%      % Test for no difference between 6-cylinder and 8-cylinder cars
%      coefTest(lme,[0 0 1 -1])
%
%   See also ANOVA, coefCI, fixedEffects, randomEffects, covarianceParameters. 

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
                        
        % (1) Save names for the big sparse Zs used to fit a standard LME.
        model.RandomInfo.ZsColNames = makeSparseZNames(model);
        
        % (2) Fit a LME model in standard form.        
        model.slme = fitStandardLMEModel(model);
       
        % (3) Fill in DFE and Coefs so that get.NumCoefficients and 
        % get.NumEstimatedCoefficients in ParametricRegression work as
        % required. Also fill in CoefficientCovariance.
        model.Coefs = model.slme.betaHat;
        N = model.NumObservations; % only contribution from subset.
        P = length(model.Coefs);
        model.DFE = N - P;
        model.CoefficientCovariance = model.slme.covbetaHat;
        
    end % end of fitter.
    
    % The predict method would normally take a dataset/table array or a matrix
    % containing all variables.  This method exists to allow prediction
    % with a matrix that contains only the required predictor variables
    % without blowing it up to contain all the variables only to then pull
    % it apart to get the design matrix.
    function ypred = predictPredictorMatrix(model,Xpred) %#ok<INUSD>
        ypred = [];
    end % end of predictPredictorMatrix.
    
    function D = get_diagnostics(model,type) %#ok<INUSD>
        D = [];
    end % end of get_diagnostics.

    function L = getlogLikelihood(model)
        
        % L will be the maximized log-likelihood if FitMethod is ML and the
        % maximized restricted log-likelihood if FitMethod is REML. We also
        % account for the effect of observation weights.
        w = model.ObservationInfo.Weights;
        subset = model.ObservationInfo.Subset;
        w = w(subset);
        L = model.slme.loglikHat + 0.5*sum(log(w));
        
    end % end of getlogLikelihood.

    function L0 = logLikelihoodNull(model)         %#ok<MANU>

        L0 = NaN;
        
    end % end of logLikelihoodNull.       
    
end

% Other inherited protected methods that we choose to redefine.
methods(Access='protected')
    
    function model = postFit(model)                
                        
        % (1) Set MSE and Response.
        model.MSE = model.slme.sigmaHat^2;        
        model.Response = response(model);
        
        % (2) Get optimized log likelihood or restricted log likelihood.
        % Also set model.LogLikelihoodNull to NaN since this quantity is
        % not well defined for a LME.
        model.LogLikelihood = getlogLikelihood(model);
        model.LogLikelihoodNull = logLikelihoodNull(model);
        
        % (3) Get SSE, SST and SSR.
        [model.SSE,model.SSR,model.SST] = getSumOfSquares(model);
       
        % (4) Get CoefficientNames.
        model.CoefficientNames = getCoefficientNames(model);                
        
        % (5) Save a table of covariance parameters in the object.
        [~,~,model.CovarianceTable] = covarianceParameters(model);
        
        % (6) Set total number of parameters.
        %model.slme.Psi.NumParametersExcludingSigma + (model.slme.p+1);
        
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
            Standardized = residuals(model,'ResidualType','standardized');

            % (2) Assemble into a table.
            r = table(Raw,Pearson,Standardized,...
                'RowNames',model.ObservationNames);        
        else
            % (1) Ensure that residualtype is a valid 'ResidualType' string.
            residualtype ...
                = LinearMixedModel.validateResidualType(residualtype);  
            
            % (2) Now get the appropriate conditional residual.
            switch lower(residualtype)                
                case 'raw'
                    r = residuals(model,'ResidualType','raw');
                case 'pearson'
                    r = residuals(model,'ResidualType','pearson');
                case 'standardized'
                    r = residuals(model,'ResidualType','standardized');                    
            end
        end
        
    end % end of get_residuals.
        
    % (7) Called by get.Fitted.
    function yfit = get_fitted(model)     
        
        % (1) Return conditional fitted values - the default in fitted.
        yfit = fitted(model);
        %yfit = predict(model);
        
    end % end of get_fitted.    
    
end

% No arg constructor.
methods(Access='public', Hidden=true)
   
    function lme = LinearMixedModel(varargin)
%   The LinearMixedModel constructor is not intended to be called directly.
%   Use FITLME to create a LinearMixedModel by fitting to data.
%
%   See also FITLME.        
        st = dbstack;
        isokcaller = false;
        if (length(st) >= 2)
            isokcaller = any(strcmpi(st(2).name,{'LinearMixedModel.fit','LinearMixedModel.fitmatrix'}));
        end
        if (nargin == 0 && isokcaller == true)
            lme.Formula = classreg.regr.LinearMixedFormula('y ~ -1');
            return;
        end        
        error(message('stats:LinearMixedModel:NoConstructor'));
    end
       
end    

% Display related private methods
methods(Access='private')        
    
    function displayHeadLine(model)
        
        % (1) Display headline.
        isLoose = strcmp(get(0,'FormatSpacing'),'loose');
        if (isLoose), fprintf('\n'); end                 
        headline = getString(message('stats:LinearMixedModel:Display_headline',model.FitMethod));
        headline = LinearMixedModel.formatBold(headline);
        fprintf('%s\n',headline);                      
        fprintf('\n');
        
    end % end of displayHeadLine.
    
    function displayFormula(model)
                
        % (1) Fit the formula string to the command window width.               
        formulaheadline = getString(message('stats:LinearMixedModel:Display_formula'));
        formulaheadline = LinearMixedModel.formatBold(formulaheadline);
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
        
        % (1) AIC, BIC etc. table.        
        modelfitstatsheadline = getString(message('stats:LinearMixedModel:Display_modelfitstats'));
        modelfitstatsheadline = LinearMixedModel.formatBold(modelfitstatsheadline);
        fprintf('%s\n',modelfitstatsheadline);        
        crittable = modelCriterionLME(model);        
        crittable = LinearMixedModel.removeTitle(crittable);
        disp(crittable);
        
    end % end of displayModelFitStats.
    
    function displayFixedStats(model)
        
        % (1) Stats table for fixed effects coefficients.        
        fixedstatsheadline = getString(message('stats:LinearMixedModel:Display_fixedstats'));
        fixedstatsheadline = LinearMixedModel.formatBold(fixedstatsheadline);
        fprintf('%s\n',fixedstatsheadline);        
        ds = model.Coefficients;        
        ds = LinearMixedModel.removeTitle(ds);
        disp(ds);
        
    end % end of displayFixedStats.
    
    function displayCovarianceStats(model)
        
        % (1) Get headline string.                  
        covariancestatsheadline = getString(message('stats:LinearMixedModel:Display_covariancestats'));
        covariancestatsheadline = LinearMixedModel.formatBold(covariancestatsheadline);
        fprintf('%s\n',covariancestatsheadline);        

        % (2) Number of grouping variables
        R = model.GroupingInfo.R;
        
        % (3) If there are random effects, extract the covariance
        % table and display it. Also, display the error covariance.
        %indent = '    ';
        lev = model.GroupingInfo.lev;
        for k = 1:(R+1)
            % (1) Name of current group.
            if k > R                
                gname = getString(message('stats:LinearMixedModel:String_error'));                
                fprintf('%s\n',[getString(message('stats:LinearMixedModel:String_group')),': ',gname]);
            else
                gname = model.GroupingInfo.GNames{k};                
                fprintf('%s\n',[getString(message('stats:LinearMixedModel:String_group')),': ',gname,' (',num2str(lev(k)),' ',getString(message('stats:LinearMixedModel:String_levels')),')']);
            end
            % (2) Delete the Group column from covtable{k}
            % if required.
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
            % (3) Display covtable{k}.
            ds = model.CovarianceTable{k};
            ds = LinearMixedModel.removeTitle(ds);
            disp(ds);
        end
        
    end % end of displayCovarianceStats.
    
    function displayModelInfo(model)
            
        % (1) Get headling string.                                      
        modelinfoheadline = getString(message('stats:LinearMixedModel:Display_modelinfo'));
        modelinfoheadline = LinearMixedModel.formatBold(modelinfoheadline);
        fprintf('%s\n',modelinfoheadline);        
        
        % (2) Number of observations.
        N = model.slme.N;
        
        % (3) Number of fixed effects.
        p = model.slme.p;
        
        % (4) Number of random effects.
        q = model.slme.q;        

        % (5) Number of covariance parameters.
        ncov = model.slme.Psi.NumParametersExcludingSigma + 1;
        
        % (6) Show the N, p, q and ncov in a dataset.
        indent = '    ';
        fprintf('%-35s %6d\n',[indent,getString(message('stats:LinearMixedModel:ModelInfo_numobs'))],N);
        fprintf('%-35s %6d\n',[indent,getString(message('stats:LinearMixedModel:ModelInfo_fecoef'))],p);
        fprintf('%-35s %6d\n',[indent,getString(message('stats:LinearMixedModel:ModelInfo_recoef'))],q);
        fprintf('%-35s %6d\n',[indent,getString(message('stats:LinearMixedModel:ModelInfo_covpar'))],ncov);
        fprintf('\n');        
        
    end % end of displayModelInfo.
    
end

% Private utility methods.
methods (Access={?classreg.regr.LinearLikeMixedModel})
                                  
    function crittable = modelCriterionLME(model)
%modelCriterionLME - Compute table containing model criterion info.
%   crittable = modelCriterionLME(model) takes a LinearMixedModel object model
%   and computes a table of model criterion such as AIC, BIC, logLik and
%   Deviance.
        
    % Create a structure stats for classreg.regr.modelutils.modelcriterion
    % such that:
    %
    % (1) stats.LogLikelihood is the maximized log likelihood or restricted 
    % log likelihood
    % 
    % (2) stats.NumCoefficients is the number of unknown coefficients in 
    % the model
    %
    % (3) stats.NumObservations is the effective number of observations 
    % used in the fit    
    
        % (1) Get N and p.
        N = model.slme.N;
        p = model.slme.p;
        
        % (2) Add the fixed effects + residual variance to NumCoefficients.
        stats.NumCoefficients = ...
            model.slme.Psi.NumParametersExcludingSigma + (p+1);
        
        % (3) Get effective number of observations based on FitMethod.
        switch lower(model.FitMethod)
            case {'ml'}
                stats.NumObservations = N;
            case {'reml'}
                stats.NumObservations = (N-p);
            otherwise                      
                error(message('stats:LinearMixedModel:BadFitMethod'));
        end
        
        % (4) Set maximized log likelihood or maximized restricted log 
        % likelihood. model.LogLikelihood also accounts for observation
        % weights.
        stats.LogLikelihood = model.LogLikelihood;
        
        % (5) Call modelcriterion utility function.
        crit = classreg.regr.modelutils.modelcriterion(stats,'all',true);
        
        % (6) Get deviance and create output table.
        Deviance = -2*stats.LogLikelihood;
        crittable = table(crit.AIC, crit.BIC, stats.LogLikelihood, Deviance,...
            'VariableNames',{'AIC' 'BIC' 'LogLikelihood' 'Deviance'});        

        % (7) Add a title to the table.        
         ttl = getString(message('stats:LinearMixedModel:Title_modelfitstats'));
         crittable = classreg.regr.lmeutils.titleddataset(crittable,ttl);
        
    end % end of modelCriterionLME.
        
    function w = getCombinedWeights(model,reduce)
%getCombinedWeights - Get the effective weights for this LinearMixedModel.        
%   w = getCombinedWeights(model,reduce) takes an object model of type
%   LinearMixedModel and returns w, the effective observation weights
%   vector. If reduce is true then the length of w equals the number of
%   observations actually used in the fit. If reduce is false then the
%   length of w equals the total number of observations in the input data
%   even if not used in the fit.        

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
%fitStandardLMEModel - Fits a LME model in standard form.
%   slme = fitStandardLMEModel(model) takes an object model of type
%   LinearMixedModel and returns a fitted StandardLinearMixedModel slme
%   such that slme is ready for stats.        

        % (1) Make big sparse Zs matrix and covariance matrix object Psi.
        Z = model.RandomInfo.Z;
        q = model.RandomInfo.q;
        lev = model.GroupingInfo.lev;
        Gid = model.GroupingInfo.Gid;
        N = model.NumObservations;
        Zs = LinearMixedModel.makeSparseZ(Z,q,lev,Gid,N);
        Psi = makeCovarianceMatrix(model);

        % (2) Exclude observations, NaNs or Infs. Already done.

        % (3) Scale input to StandardLinearMixedModel using obs weights.
        reduce = true;
        
        w = getCombinedWeights(model,reduce);
        Xw =scaleMatrixUsingWeights(model.FixedInfo.X,w);
        yw = scaleMatrixUsingWeights(model.y,w);
        Zsw =scaleMatrixUsingWeights(Zs,w);

        % (4) Fit a StandardLinearMixedModel and make it ready for stats.
        dofit = true;
        dostats = false;
        
        
        slme = classreg.regr.lmeutils.StandardLinearMixedModel(Xw,yw,Zsw,Psi,model.FitMethod,...
            dofit,dostats,'Optimizer',model.Optimizer,...
            'OptimizerOptions',model.OptimizerOptions,...
            'InitializationMethod',model.StartMethod,...
            'CheckHessian',model.CheckHessian);
        
        B0 = 1e6*ones(size(slme.betaHat));
        iter=1; maxiter=20;
        while norm(slme.betaHat-B0)/norm(B0) > 1e-2 && iter < maxiter
            
            B0 =slme.betaHat;
            
            resid = yw - Xw*slme.betaHat - Zsw*slme.thetaHat;
            a = nirs.math.ar_fit(resid, model.pmax);
            
            % create a whitening filter from the coefficients
            f = [1; -a(2:end)];
            
            % filter the design matrix
            Xf = myFilter(f,Xw);
            Zf = myFilter(f,Zsw);
            yf = myFilter(f,yw);
            
            for iter=1:10
                resid = yf - Xf*slme.betaHat - Zf*slme.thetaHat;
                w=wfun(resid);
                
                Xfw=diag(w)*Xf;
                yfw=diag(w)*yf;
                Zfw=diag(w)*Zf;
                slme = classreg.regr.lmeutils.StandardLinearMixedModel(Xfw,yfw,Zfw,Psi,model.FitMethod,...
                    dofit,dostats,'Optimizer',model.Optimizer,...
                    'OptimizerOptions',model.OptimizerOptions,...
                    'InitializationMethod',model.StartMethod,...
                    'CheckHessian',model.CheckHessian);
                
            end
            dostats = true;
            slme = classreg.regr.lmeutils.StandardLinearMixedModel(Xfw,yfw,Zfw,Psi,model.FitMethod,...
                dofit,dostats,'Optimizer',model.Optimizer,...
                'OptimizerOptions',model.OptimizerOptions,...
                'InitializationMethod',model.StartMethod,...
                'CheckHessian',model.CheckHessian);
            
            
        end
         fprintf(1,'.');
        
        
    end % end of fitStandardLMEModel.
            
    function [SSE,SSR,SST] = getSumOfSquares(model)
%getSumOfSquares - Computes SSE, SSR and SST.
%   [SSE,SSR,SST] = getSumOfSquares(model) takes an object model of type
%   LinearMixedModel and computes SSE, SSR and SST using conditional fitted
%   values.
%
%   If Y is the response vector and F is a vector of conditional fitted 
%   values then:
%
%   SSE = sum( (Y-F).^2 );
%   SSR = sum( (F-mean(F)).^2 );
%   SST = sum( (Y-mean(Y)).^2 );
%
%   Only observations actually used in the fit should be used for
%   calculating SSE, SSR and SST. If observation weights have been supplied
%   then: 
%
%   SSE = sum( w.* (Y-F).^2 );
%   SSR = sum( w.* (F-meanw(F)).^2 );
%   SST = sum( w.* (Y-meanw(Y)).^2 );
%
%   where meanw(Y) is the weighted mean of Y using weights w. Same for
%   meanw(F). If the model has an intercept then meanw(Y) should be equal
%   to meanw(F).

        % (1) Get *conditional* fitted values.
        F = fitted(model);
        
        % (2) Get response vector.
        Y = response(model);
        
        % (3) Get the weight vector.
        w = model.ObservationInfo.Weights;
        
        % (4) Subset of observations actually used in the fit.
        subset = model.ObservationInfo.Subset;
        
        % (5) Keep only observations used in the fit in F, Y and w.
        F = F(subset);
        Y = Y(subset);
        w = w(subset);        
        
        % (6) Weighted mean of Y and F.
        Y_mean_w = sum(w.*Y)/sum(w);
        %F_mean_w = sum(w.*F)/sum(w);
        
        % (7) Compute SSE.
        SSE = sum(w.*((Y - F).^2));   
        
        % (8) Compute SSR. Y_mean_w should be the same as F_mean_w if the
        % model has an intercept.
        SSR = sum(w.*((F - Y_mean_w).^2));
        
        % (9) Compute SST.
        SST = SSE + SSR;
        %SST = sum(w.*((Y - Y_mean_w).^2));                            

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
%   plus 1 (for residual noise variance).
       
        % (1) Get number of fixed effects parameters.
        numfixpar = model.slme.p;
        
        % (2) Get number of covariance parameters excluding noise std.
        numcovpar = model.slme.Psi.NumParametersExcludingSigma;
        
        % (3) Get np.
        np = numfixpar + numcovpar + 1;
        
    end % end of getTotalNumberOfParameters.
    
    function ysim = randomForSLRT(model,S)
%randomForSLRT - Helper method for generating random data in SLRT.
%   ysim = randomForSLRT(model,S) takes a LinearMixedModel model and a
%   RandStream object S and generates random data from the fitted model at
%   the original design points. The effect of originally supplied
%   observation weights is also included. The length of ysim is equal to
%   the length of model.y or model.slme.y. In other words, excluded
%   observations are not included in ysim.

        % (1) Get the slme object from model.
        slmeObj = model.slme;

        % (2) First get W^(1/2)*ysim and then get ysim from that. Pass
        % along the RandStream object S to the random method of
        % StandardLinearMixedModel.
        ysim = random(slmeObj,S,slmeObj.X,slmeObj.Z);
        reduce = true;
        w = getCombinedWeights(model,reduce);
        ysim = ysim ./ sqrt(w);                       
        
    end % end of randomForSLRT.
    
    function loglik = doFitFunForSLRT(model,ysim)
%doFitFunForSLRT - Helper method for fitting the model to new data in SLRT.
%   loglik = doFitFunForSLRT(model,ysim) takes a fitted LinearMixedModel
%   model and a new realization ysim of model.y and fits the model to ysim.
%   The maximized log likelihood or the maximized restricted log likelihood
%   is returned. ysim is on the scale of the original data i.e., does not
%   include the effect of observation weights.
        
        % (1) Prepare response to be passed on to model.slme. We need to
        % scale ysim by the square root of the weights before passing it on
        % to model.slme.
        reduce = true;
        w = getCombinedWeights(model,reduce);
        yw = LinearMixedModel.scaleMatrixUsingWeights(ysim,w);

        % (2) Get the slme object from model.
        slmeObj = model.slme;
        
        % (3) Modify the slme object by setting its y to yw.
        slmeObj.y = yw; 
        
        % (4) Refit the model.
        slmeObj = refit(slmeObj); 
        
        % (5) Get the maximized log likelihood or the maximized restricted 
        % log likelihood. Include the contribution of observation weights.
        loglik = slmeObj.loglikHat + 0.5*sum(log(w));         
        
    end % end of doFitFunForSLRT.
        
end


% Private static utility methods.
methods(Static,Access='private')
       
    function [X,XColNames,XCols2Terms] = designFromFormula(F,ds,dummyvarcoding)
%designFromFormula - Compute a design matrix from a LinearFormula object.
%   [X,XColNames,XCols2Terms] = designFromFormula(F,ds,dummyvarcoding)
%   takes a LinearFormula object F, a dataset/table array ds, and a string
%   dummyvarcoding specifying the dummy variable coding to use and returns
%   the N-by-P design matrix X, a cell array XColNames of length P
%   containing the name of each column in X and an integer array
%   XCols2Terms of size 1-by-P indicating which "term" a particular column
%   in X belongs to. All columns of X that belong to the same "term" can be
%   entered into an ANOVA table with an F-test.
        
        [X,~,~,XCols2Terms,XColNames] = ...
            classreg.regr.modelutils.designmatrix(ds,'Model',F.Terms,...
            'PredictorVars',F.PredictorNames,'ResponseVar',F.ResponseName,...
            'DummyVarCoding',dummyvarcoding);
        
        
%         [design,~,~,coefTerm,coefNames] ...
%                     = classreg.regr.modelutils.designmatrix(X,'Model',terms(:,varLocs), ...
%                     'DummyVarCoding',model.DummyVarCoding, ...
%                     'CategoricalVars',model.VariableInfo.IsCategorical(varLocs), ...
%                     'CategoricalLevels',model.VariableInfo.Range(varLocs));
        
        
        
    end % end of designFromFormula.
    
 function Xw = scaleMatrixUsingWeights(X,w)
%scaleMatrixUsingWeights - Multiplies each row of X using sqrt(w).
%   Xw = scaleMatrixUsingWeights(X,w) takes a N-by-P matrix X and a N-by-1 
%   weight vector w and outputs a scaled matrix Xw such that: 
%           Xw(i,:) = X(i,:) * sqrt(w(i)).       
        
        % (1) Ensure that X is a numeric, real matrix.
        assert( isnumeric(X) & isreal(X) & ismatrix(X) );
        
        % (2) Get the size of X.
        N = size(X,1);
        
        % (3) Ensure that w is a weight vector of length N.
        w = LinearMixedModel.validateWeights(w,N);

        % (4) Form Xw. If all weights are 1 just copy X into Xw.
        if all(w == 1)
            Xw = X;
        else
            Xw = bsxfun(@times,X,sqrt(w));
        end

    end % end of scaleMatrixUsingWeights.
          
    function [ds,formula] = convertXYZGToDataset(X,Y,Z,G,XNames,YName,ZNames,GNames)
%convertXYZGToDataset - convert matrix to table representation.      
%   [ds,formula] = convertXYZGToDataset(X,Y,Z,G,XNames,YName,ZNames,GNames)        
%   takes validated inputs X, Y, Z, G, XNames, YName, ZNames and GNames and
%   creates a table and formula representation of the LME model.
%
%   X - N-by-P matrix
%   Y - N-by-1 vector
%   Z - R-by-1 cell array. Z{k} is a N-by-q(k) matrix.
%   G - R-by-1 cell array. G{k} is a length N grouping variable.
%   XNames - P-by-1 cell array of strings containing names of columns of X.
%    YName - String containing name for Y.
%   ZNames - R-by-1 cell array. ZNames{k} is a q(k)-by-1 cell array of 
%            strings containing names for columns of Z{k}.
%   GNames - R-by-1 cell array of strings. GNames{k} contains the name of
%            grouping variable G{k}.

        % (1) Assume that inputs are validated.
        
        % (2) Create empty table.
        ds = table();
        
        % (3) First put columns of X into ds.
        ds = LinearMixedModel.addColumnsToDataset(ds,X,XNames);
       
        % (4) Put Y in ds.
        ds = LinearMixedModel.addColumnsToDataset(ds,Y,{YName});
        
        % (5) Loop over cell array Z and put each element in ds.
        R = length(Z);
        for k = 1:R
            ds = LinearMixedModel.addColumnsToDataset(ds,Z{k},ZNames{k});            
        end
        
        % (6) Loop over elements of G and put each element in ds.
        for k = 1:length(G)
            if ~ismember(GNames{k},ds.Properties.VariableNames)
                ds.(GNames{k}) = G{k};
            end
        end
        
        % (7) Create a formula string with no predictors.
        formula = YName;

        % (8) Add fixed effects spec to formula.
        fespec = LinearMixedModel.getFixedRandomSpec(XNames,[],false);
        formula = [formula,' ~ ',fespec];

        % (9) Create R random effects specs and add to formula.
        for k = 1:R
            respec = LinearMixedModel.getFixedRandomSpec(ZNames{k},...
                GNames{k},false);
            formula = [formula,' + ',respec];             %#ok<AGROW>
        end        
        
    end % end of convertXYZGToDataset.
    
    function ds = addColumnsToDataset(ds,X,XNames)
%addColumnsToDataset - Add columns of X to ds.
%   ds = addColumnsToDataset(ds,X,XNames) takes a table array ds, a N by
%   P matrix X, a length P cell array of strings XNames and adds the
%   columns of X to ds. Only those columns will be added to ds which are
%   not already present in ds.
        
        % (1) Assume X and XNames have been validated.
        
        % (2) Add columns of X to ds.
        p = size(X,2);
        for k = 1:p
            % We are trying to add a column with name XNames{k} to ds. Do
            % this only if ds does not already have this column.
            if ~ismember(XNames{k},ds.Properties.VariableNames)
                ds.(XNames{k}) = X(:,k);
            end
        end

    end % end of addColumnsToDataset.        
        
    function spec = getFixedRandomSpec(varnames,gname,wantintercept)
%getFixedRandomSpec -Get fixed or random effect specification.
%   spec = getFixedRandomSpec(varnames,gname) takes a cell array of strings
%   varnames and a string gname and returns a random effects specification
%   corresponding to variables in varnames and grouping variable in gname.
%   If gname is [], a fixed effects specification corresponding to
%   variables in varnames is returned. wantintercept is either true or
%   false. If wantintercept is false then a '-1' is also added to the spec.
%   Default value of gname is [] and default value of wantintercept =
%   false.
%
%   Example: varnames = {'x1','x2'}, gname = 'g'. 
%       if wantintercept = true : spec = '(1 + x1 + x2 | g)'
%       if wantintercept = false: spec = '(-1 + x1 + x2 | g)'
%
%   Example: varnames = {'x1','x2','x3'}, gname = [].
%       if wantintercept = true : spec = '(1 + x1 + x2 + x3)'
%       if wantintercept = false: spec = '(-1 + x1 + x2 + x3)'
        
        % (1) Set default values.
        switch nargin
            case 1
                % gname and wantintercept not given.
                gname = [];
                wantintercept = false;
            case 2
                % wantintercept not given.
                wantintercept = false;
        end            
        
        % (2) Initialize spec.
        if wantintercept == true
            spec = '1';
        else
            spec = '-1';            
        end
        
        % (3) Fixed effects specification.
        p = length(varnames);        
        for k = 1:p
            spec = [spec,' + ',varnames{k}]; %#ok<AGROW>
        end
        
        % (4) If gname is given, create a random effects spec.
        if ~isempty(gname)
            spec = ['(',spec,' | ',gname,')'];             
        end                

    end % end of getFixedRandomSpec.        
                                   
    function [lrt,siminfo] = simulatedLRT(smallModel,bigModel,smallModelName,bigModelName,nsim,alpha,options)
%simulatedLRT - Simulated likelihood ratio test.
%   [table,siminfo] = simulatedLRT(smallModel,bigModel,smallModelName,bigModelName,nsim,alpha,options)
%   takes smallModel and bigModel - two LinearMixedModel objects.
%   smallModel has name smallModelName and bigModel has name bigModelName.
%   We do a simulated likelihood ratio test under the assumption that
%   observed data arises from the simpler model smallModel. We assume that
%   nesting requirements have already been checked. nsim is the number of
%   simulations to run. alpha is used to decide the coverage of confidence
%   interval for the estimated p-value. Example, alpha = 0.05 for 95%
%   confidence interval. options is a structure containing parallel
%   processing options. Assume that nsim, alpha and options have already
%   been validated.
%
%   The output lrt should be a dataset that looks like this:
%
% Model     DF   AIC     BIC     LogLik    LRStat   pValue   Lower   Upper
% lme       4    149.22  156.17  -70.609
% altlme    6    149.43  159.85  -68.714   3.7896  0.15035     x       x


        % (1) Get useParallel and RNGscheme to be used by smarForSliceout.
        [useParallel,RNGscheme] ...
            = internal.stats.parallel.processParallelAndStreamOptions(options,true);
        
        % (2) Turn display off on smallModel.slme and bigModel.slme before
        % starting simulations.
        smallModel.slme = turnOffOptimizerDisplay(smallModel.slme);
          bigModel.slme = turnOffOptimizerDisplay(bigModel.slme);
        
        % (3) Make the loopbody function.
        loopbody = LinearMixedModel.makeloopbodyFunSLRT(smallModel,bigModel);
        
        % (4) Call smartForSliceout.
        TH0 ...
            = internal.stats.parallel.smartForSliceout(nsim,loopbody,useParallel,RNGscheme);

        % (5) Observed LR test statistic.
        T = 2*(bigModel.LogLikelihood - smallModel.LogLikelihood);

        % (6) Compute the simulation based p-value and its 100*(1-alpha)
        % confidence interval:
        %       pvalueSim = ( 1 + sum( TH0 >= T ) ) / ( 1 + nsim ).
        [pvalueSim,pvalueSimCI] ...
            = binofit( 1+sum(TH0 >= T), 1+nsim, alpha );
        
        % (7) Initialize lrt to be the standard LRT table.
        lrt = LinearMixedModel.standardLRT(smallModel,bigModel,...
            smallModelName,bigModelName);               
        
        % (8) At this point, lrt looks like this:
        % Model     DF   AIC     BIC     LogLik    LRStat  deltaDF  pValue
        % lme       4    149.22  156.17  -70.609
        % altlme    6    149.43  159.85  -68.714   3.7896    2      0.15035
        % Delete the deltaDF column from lrt.
        lrt.deltaDF = [];
        
        % (9) Modify the pValue column so that it contains the simulated
        % p-value.
        pValue = zeros(2,1);
        pValue(1) = 0;
        pValue(2) = pvalueSim;
        pValueAbsent = [true;false];
        lrt.pValue = internal.stats.DoubleTableColumn(pValue,pValueAbsent);
        
        % (10) Set confidence interval for simulated p-value.
            Lower = zeros(2,1);
            Lower(1) = 0;
            Lower(2) = pvalueSimCI(1);
            LowerAbsent = [true;false];
            lrt.Lower ...
                = internal.stats.DoubleTableColumn(Lower,LowerAbsent);

            Upper = zeros(2,1);
            Upper(1) = 0;
            Upper(2) = pvalueSimCI(2);
            UpperAbsent = [true;false];
            lrt.Upper ...
                = internal.stats.DoubleTableColumn(Upper,UpperAbsent);
        
        % (11) Add a title to lrt.        
        ttl = getString(message('stats:LinearMixedModel:Title_SLRT',num2str(nsim),num2str(alpha)));
        lrt = classreg.regr.lmeutils.titleddataset(lrt,ttl);  
            
        % (12) Save additional simulation related info in siminfo.
            % (1) Number of simulations.
            siminfo.nsim = nsim;
            % (2) Selected value for alpha.
            siminfo.alpha = alpha;
            % (3) Simulated p-value.
            siminfo.pvalueSim = pvalueSim;
            % (4) 100*(1-alpha) confidence interval on pvalueSim.
            siminfo.pvalueSimCI = pvalueSimCI;
            % (5) Number of free parameters in bigModel minus those in
            % smallModel.
                DF = lrt.DF;
                siminfo.deltaDF = DF(2) - DF(1);
            % (6) Simulated distribution of T under the null hypothesis
            % that smallModel adequately explains the observed data.
            siminfo.TH0 = TH0;
      
    end % end of simulatedLRT.
    
    function fun = makeloopbodyFunSLRT(smallModel,bigModel)
        
        % Helper function to make loop body function for smartForSliceout
        % in SLRT.
        fun = @loopbodyFun;
        
        function TH0 = loopbodyFun(~, S)

            % (1) Get ysim from smallModel. S is a RandStream object.
            ysim = randomForSLRT(smallModel,S);
            
            % (2) Get LR statistic for ysim.
            loglik1 = doFitFunForSLRT(smallModel,ysim);            
            loglik2 = doFitFunForSLRT(bigModel,ysim);            
            TH0 = 2*(loglik2 - loglik1);
            
        end % end of loopbodyFun.        
        
    end
            
end

methods(Static, Access='protected')
    
    function checkNestingRequirement(smallModel,bigModel,smallModelName,bigModelName,isSimulatedTest)
%checkNestingRequirement - Checks if smallModel is nested in bigModel.       
%  checkNestingRequirement(smallModel,bigModel,smallModelName,bigModelName,isSimulatedTest)
%  takes two LinearMixedModel objects, smallModel and bigModel. smallModel
%  is the potentially smaller model that is nested in the bigger model
%  bigModel. smallModel has name smallModelName and bigModel has name
%  bigModelName. We try to check if smallModel is really nested in bigModel
%  or not. This technique is not fool proof since we do not check the
%  covariance structure of random effects but it is better than no check.
%  isSimulatedTest is either true or false. If isSimulatedTest is true then
%  the intention is to check nesting for a simulated likelihood ratio test.
%  If isSimulatedTest is false then the intention is to check nesting for a
%  standard likelihood ratio test.
%
%   What is checked?
%
%   If isSimulatedTest is true (simulated LRT)
%   
%  (1) smallModel and bigModel have been fit using the same FitMethod.
%
%  (2) smallModel and bigModel have the same response vector (used in fit).
%
%   Suppose Xsmall is the fixed effects design matrix of smallModel in
%   standard form and Xbig is the fixed effects design matrix of bigModel
%   in standard form actually used in the fit. Xsmall and Xbig are not 
%   scaled by weights. 
%
%  (3) If FitMethod is REML, span of Xbig and Xsmall should be the same.
%
%  (4) The maximized log likelihood (or restricted log likelihood) of
%  bigModel must be >= the maximized log likelihood of the smallModel.
%
%   If isSimulatedTest is false (standard LRT) then in addition to (1), (2)
%   (3) and (4), we also check:
%
%  (5) smallModel and bigModel have the same weights vector (used in fit).
%
%  (6) If FitMethod is ML, the span of Xbig must contain Xsmall.
%
%   Suppose Zsmall is the overall random effects design matrix of
%   smallModel in standard form and Zbig is the overall random effects
%   design matrix of bigModel in standard form actually used in the fit.
%   Since weights are the same for standard LRT, we can use either the
%   weighted or unweighted Zsmall and Zbig. For convenience, we use the
%   weighted versions.
%
%  (7) The span of Zbig must contain Zsmall.

        % (1) Ensure that smallModel and bigModel are LinearMixedModel
        % objects.
        assert( isa(smallModel,'LinearMixedModel') );
        assert( isa(  bigModel,'LinearMixedModel') );
        
        % (2) Ensure that smallModelName and bigModelName are strings.
        assert( internal.stats.isString(smallModelName) );
        assert( internal.stats.isString(  bigModelName) );

        % (3) Ensure that isSimulatedTest is a scalar logical.
        assert( isscalar(isSimulatedTest) & islogical(isSimulatedTest) );
                
        % (4) Have smallModel and bigModel been fit using the same
        % FitMethod?
        fitmethodsmall = smallModel.FitMethod;
        fitmethodbig   =   bigModel.FitMethod;                
        assertThat(isequal(fitmethodsmall,fitmethodbig),'stats:LinearMixedModel:NestingCheck_fitmethod',smallModelName,bigModelName);
                
        % (5) Ensure that smallModel and bigModel have the same response 
        % vector, actually used in the fit.       
        ysmall = smallModel.y;
        ybig   =   bigModel.y;                  
        assertThat(isequaln(ysmall,ybig),'stats:LinearMixedModel:NestingCheck_response',smallModelName,bigModelName);
        
        % (6) Get Xbig and Xsmall.
        Xsmall = smallModel.FixedInfo.X;
        Xbig   =   bigModel.FixedInfo.X;
        
        % (7) The span of Xbig and Xsmall must be identical if FitMethod is
        % REML.
        if strcmpi(smallModel.FitMethod,'reml')            
            assertThat(isequaln(Xsmall,Xbig),'stats:LinearMixedModel:NestingCheck_spanX',smallModelName,bigModelName);            
        end
        
        % (8) Is the maximized log likelihood (or restricted log
        % likelihood) of bigModel >= that of smallModel?
        logliksmall = smallModel.LogLikelihood;
        loglikbig   =   bigModel.LogLikelihood;               
        assertThat(loglikbig >= logliksmall,'stats:LinearMixedModel:NestingCheck_loglik',bigModelName,smallModelName);
        
        if isSimulatedTest == false
           % Additional checks. 
            
           % (9) Ensure the smallModel and bigModel have the same weights
           % vector used in the fit.
           wsmall = getCombinedWeights(smallModel,true);
           wbig   = getCombinedWeights(  bigModel,true);                     
           assertThat(isequaln(wsmall,wbig),'stats:LinearMixedModel:NestingCheck_weights',smallModelName,bigModelName);
        
           if strcmpi(smallModel.FitMethod,'ml')           
               % (10) The span of Xbig must contain Xsmall.                               
               assertThat(LinearMixedModel.isMatrixNested(Xsmall,Xbig),'stats:LinearMixedModel:NestingCheck_nestedspanX',smallModelName,bigModelName);               
           end
               
           % (11) The span of Zbig must contain Zsmall.
           Zsmall = smallModel.slme.Z;
           Zbig   =   bigModel.slme.Z;                      
           assertThat(LinearMixedModel.isMatrixNested(Zsmall,Zbig),'stats:LinearMixedModel:NestingCheck_nestedspanZ',smallModelName,bigModelName);                      
        end
        
    end % end of checkNestingRequirement.
    
end

% Private, static utility methods for input validation.
methods(Static, Access='protected')    
               
    function fitmethod = validateFitMethod(fitmethod)
%validateFitMethod - Validate the 'FitMethod' parameter.
%   fitmethod = validateFitMethod(fitmethod) takes a potential 'FitMethod'
%   parameter and either returns the validated parameter fitmethod or 
%   throws an error message. 
%
%   What is checked?
%
%   (1) fitmethod is a string containing either 'ML' or 'REML'.
%
%   If (1) is *not* satisfied, then an error message is thrown.
        
        fitmethod = internal.stats.getParamVal(fitmethod,...
            LinearMixedModel.AllowedFitMethods,'FitMethod');

    end % end of validateFitMethod.
    
    function w = validateWeights(w,N)
%validateWeights - Validates the 'Weights' parameter.        
%   w = validateWeights(w,N) takes a potential 'Weights' parameter w and an
%   integer N specifying the expected size of w. Either the output w is a
%   validated column vector of weights or an error is thrown.
%
%   What is checked?
%
%   (1) w is a numeric, real, vector of length N.
%
%   (2) All elements of w are not-NaN, > 0 and < Inf.
%
%   If any of (1) and (2) do not hold, then an error message is thrown. 

        % (1) N must be >= 0 and it must be an integer.
        assert(N >= 0 & internal.stats.isScalarInt(N));        
        
        % (2) (a) w must be a numeric, real, vector with length(w) == N. 
        %     (b) All elements of w are non-NaN, >= 0 and < Inf.        
        assertThat(isnumeric(w) & isreal(w) & isvector(w) & length(w)==N,'stats:LinearMixedModel:BadWeights',num2str(N));
        assertThat(      all( ~isnan(w) & w >= 0 & w < Inf )            ,'stats:LinearMixedModel:BadWeights',num2str(N));        
   
        % (3) Make w into a column vector if required.
        if size(w,1) == 1
            w = w';
        end
        
    end % end of validateWeights.
            
    function options = validateOptions(options)
%validateOptions - Validate the 'Options' parameter.         
%   options = validateOptions(options) accepts a tentative 'Options' 
%   parameter and returns the validated options parameter or throws an 
%   error message.
%
%   What is checked?
%
%   (1) options must be a structure.
%
%   If (1) is *not* satisfied, an error message is thrown.
               
        % (1) Check if options is a structure. If not, throw an error.        
        assertThat(isstruct(options),'stats:LinearMixedModel:MustBeStruct');

    end % end of validateOptions.
    
    function [optimizer,optimizeroptions] = ...
            validateOptimizerAndOptions(optimizer,optimizeroptions)
%validateOptimizerAndOptions - Validate the 'Optimizer' parameter and the 'OptimizerOptions' parameter.
%   [optimizer,optimizeroptions] = validateOptimizerAndOptions(optimizer,optimizeroptions)
%   takes a potential string optimizer and a structure or an object
%   optimizeroptions for the optimizer and returns the validated values.
%
%   What is checked?
%
%   (1) optimizer is a string that occurs in the cell array
%   LinearMixedModel.AllowedOptimizers.
%
%   (2) If optimizer is 'quasinewton', optimzeroptions can either be empty
%   or a struct.
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

        % (1) Ensure that optimizer is sensible.
        optimizer = internal.stats.getParamVal(optimizer,...
            LinearMixedModel.AllowedOptimizers,'Optimizer');        
        
        % (2) Check for license if required.
        switch lower(optimizer)
            case 'quasinewton'                                                            
                % (1) No license check needed.      
                
            case {'fminunc'}                
                % (1) Check for license when using fminunc.
                if ( license('test','Optimization_Toolbox') == false )                    
                    error(message('stats:LinearMixedModel:LicenseCheck_fminunc'));                    
                end
        end
        
        % (3) Do the following:
        %   (a) Ensure that optimizeroptions is either empty or a struct or 
        %   object created using statset or optimoptions depending on the 
        %   string optimizer and 
        %
        %   (b) Create a completely filled out optimizeroptions structure 
        %   or object as required. Either use the defaults or merge the 
        %   supplied options with the defaults.
        switch lower(optimizer)
            case 'quasinewton'                                                                            
                % (1) optimizeroptions can be empty or a struct, otherwise 
                % we have an error.                                
                assertThat(isempty(optimizeroptions) || isstruct(optimizeroptions),'stats:LinearMixedModel:OptimizerOptions_qn'); 
                
                % (2) Filled optimizeroptions.
                dflts = statset('linearmixedmodel');
                if isempty(optimizeroptions)
                    optimizeroptions = dflts;
                else
                    optimizeroptions = statset(dflts,optimizeroptions);
                end
                
            case {'fminunc'}                               
                % (1) optimizeroptions can be empty or of class
                % optim.options.SolverOptions, otherwise we have an error.                                
                assertThat(isempty(optimizeroptions) || isa(optimizeroptions,'optim.options.SolverOptions'),'stats:LinearMixedModel:OptimizerOptions_fminunc');        
                
                % (2) Filled optimizeroptions.
                if isempty(optimizeroptions)
                    optimizeroptions = optimoptions('fminunc');
                    optimizeroptions.Algorithm = 'quasi-newton';
                else
                    optimizeroptions = optimoptions('fminunc',optimizeroptions);
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
%   (1) startmethod is a string that occurs in the cell array
%   LinearMixedModel.AllowedStartMethods.
%
%   If (1) is *not* satisfied, then an error message is thrown.

        startmethod = ...
            internal.stats.getParamVal(startmethod,LinearMixedModel.AllowedStartMethods,'StartMethod');     
            
    end % end of validateStartMethod.
    
    function dfmethod = validateDFMethod(dfmethod)
%validateDFMethod - Validates the degrees of freedom method.        
%   dfmethod = validateDFMethod(dfmethod) accepts a potential degrees of 
%   freedom and validates it. If not valid, an error message is thrown.
%
%   What is checked?
%
%   (1) dfmethod is a string contained in LinearMixedModel.AllowedDFMethods
%
%   If (1) is not true, an error message is thrown.

        dfmethod = internal.stats.getParamVal(dfmethod,...
            LinearMixedModel.AllowedDFMethods,'DFMethod');     
        
    end % end of validateDFMethod.
                  
    function residualtype = validateResidualType(residualtype)
%validateResidualType - Validates the 'ResidualType' parameter.
%   residualtype = validateResidualType(residualtype) takes a potential
%   value of 'ResidualType' parameter and validates it. If not valid, an
%   error is thrown.
%
%   What is checked?
%
%   (1) residualtype must be a string contained in 
%   LinearMixedModel.AllowedResidualTypes
%
%   If (1) is not satisfied, an error message is thrown.

        residualtype = internal.stats.getParamVal(residualtype,...
            LinearMixedModel.AllowedResidualTypes,'ResidualType');
        
    end % end of validateResidualType.
            
    function prediction = validatePrediction(prediction)
%validatePrediction - Validates the 'Prediction' input.
%   prediction = validatePrediction(prediction) accepts a potential value
%   for 'Prediction', validates it and returns it. If not valid, an error
%   message is thrown.
%
%   What is checked?
%
%   (1) prediction must be a string and must be either 'curve' or 
%   'observation'.
%
%   If (1) is not satisfied, then an error is thrown.
     
        prediction = internal.stats.getParamVal(prediction,...
            {'curve','observation'},'Prediction');
        
    end % end of validatePrediction.         
                                
    function [X,Z,G] = validateXZG(X,Z,G,XName,ZName,GName)
%validateXZG - Validates inputs X, Z and G.
%   [X,Z,G] = validateXZG(X,Z,G,XName,ZName,GName) takes inputs X, Z, G
%   accepted by methods that accept a matrix form input and returns the
%   validated X, Z and G. XName, ZName and GName are the names for X, Z and
%   G respectively to be used in displaying the error messages.
%
%   What is checked?
%
%   (1) X must be a numeric, real matrix.
%
%   (2) Z must be a numeric, real matrix with size(X,1) rows or Z must be a
%   R-by-1 cell array of numeric, real matrices each with size(X,1) rows.
%
%   (3) G must be a grouping variable (cell, char, logical, categorical or
%   a numeric, real vector) of length size(X,1) or a R-by-1 cell array such
%   that each element of G is a grouping variable of length size(X,1).        
%
%   The returned X, Z and G will have the following dimensions:
%
%   1. X will be a N-by-P matrix.
%   2. Z will be a R-by-1 cell array with Z{k} a N-by-Q(k) matrix.
%   3. G will be a R-by-1 cell array with G{k} a N rows grouping variable.

        % (1) Ensure that X is sensible. Get the number of rows in X.
        X = LinearMixedModel.validateMatrix(X,XName);
        N = size(X,1);                
        
        % (2) Default value for Z.
        if isempty(Z)
            Z = zeros(N,0);
        end
            
        % (3) Ensure that Z is sensible. If Z is not a cell array, make it
        % into a column cell array. Then ensure Z{k}'s have N rows.
        if ~iscell(Z)
            Z = {Z};
        end       
        if size(Z,1) == 1
            Z = Z';
        end
        Z = LinearMixedModel.validateCellVector(Z,'Z');
        R = size(Z,1);
        for k = 1:R
            ZNamek = [ZName,'{',num2str(k),'}'];
            Z{k} = LinearMixedModel.validateMatrix(Z{k},ZNamek,N); 
        end
         
        % (4) Default value for G - R-by-1 cell array with each element
        % equal to ones(N,1).
        if isempty(G)
            G = cell(R,1);
            G(1:R) = {ones(N,1)};
        end
            
        % (5) Ensure that G is sensible. If G is not a cell array, make it
        % into a column cell array. Then validate each element. Ensure that
        % the length of G cell array is the same as that of Z cell array.
        if ~iscell(G)
            G = {G};
        end
        if size(G,1) == 1
            G = G';
        end
        G = LinearMixedModel.validateCellVector(G,'G');               
        assertThat(length(G) == R,'stats:LinearMixedModel:MustBeCellArraysOfSameLength',ZName,GName);
        for k = 1:R
            GNamek = [GName,'{',num2str(k),'}'];
            G{k} = LinearMixedModel.validateGroupingVar(G{k},GNamek,N);
        end
        
    end % validateXZG.
    
    function [X,Y,Z,G] = validateXYZG(X,Y,Z,G,XName,YName,ZName,GName)
%validateXYZG - Validates inputs X, Y, Z and G.
%   [X,Y,Z,G] = validateXYZG(X,Y,Z,G,XName,YName,ZName,GName) takes inputs
%   X, Y, Z, G accepted by fitmatrix and other methods that accept a matrix
%   form input and returns the validated X, Y, Z and G. XName, YName, ZName
%   and GName are the names for X, Y, Z and G respectively to be used in
%   displaying the error messages.
%
%   What is checked?
%
%   (1) X must be a numeric, real matrix.
%
%   (2) Y must be a numeric, real vector of length size(X,1).
%
%   (3) Z must be a numeric, real matrix with size(X,1) rows or Z must be a
%   R-by-1 cell array of numeric, real matrices each with size(X,1) rows.
%
%   (4) G must be a grouping variable (cell, char, logical, categorical or
%   a numeric, real vector) of length size(X,1) or a R-by-1 cell array such
%   that each element of G is a grouping variable of length size(X,1).        
%
%   The returned X, Y, Z and G will have the following dimensions:
%
%   1. X will be a N-by-P matrix.
%   2. Y will be a N-by-1 vector.
%   3. Z will be a R-by-1 cell array with Z{k} a N-by-Q(k) matrix.
%   4. G will be a R-by-1 cell array with G{k} a N rows grouping variable.

        % (1) Ensure that X is sensible. Get the number of rows in X.
        X = LinearMixedModel.validateMatrix(X,XName);
        N = size(X,1);
        
        % (2) Ensure that Y is sensible. First ensure that Y is a numeric,
        % real matrix. Try to make it into a column vector and then ensure 
        % that it is of size N-by-1.
        Y = LinearMixedModel.validateMatrix(Y,YName);
        if size(Y,1) == 1
            Y = Y';
        end
        Y = LinearMixedModel.validateMatrix(Y,YName,N,1);  
        
        % (3) Default value for Z.
        if isempty(Z)
            Z = zeros(N,0);
        end
            
        % (4) Ensure that Z is sensible. If Z is not a cell array, make it
        % into a column cell array. Then ensure Z{k}'s have N rows.
        if ~iscell(Z)
            Z = {Z};
        end       
        if size(Z,1) == 1
            Z = Z';
        end
        Z = LinearMixedModel.validateCellVector(Z,'Z');
        R = size(Z,1);
        for k = 1:R
            ZNamek = [ZName,'{',num2str(k),'}'];
            Z{k} = LinearMixedModel.validateMatrix(Z{k},ZNamek,N); 
        end
         
        % (5) Default value for G - R-by-1 cell array with each element
        % equal to ones(N,1).
        if isempty(G)
            G = cell(R,1);
            G(1:R) = {ones(N,1)};
        end
            
        % (6) Ensure that G is sensible. If G is not a cell array, make it
        % into a column cell array. Then validate each element. Ensure that
        % the length of G cell array is the same as that of Z cell array.
        if ~iscell(G)
            G = {G};
        end
        if size(G,1) == 1
            G = G';
        end
        G = LinearMixedModel.validateCellVector(G,'G');                
        assertThat(length(G) == R,'stats:LinearMixedModel:MustBeCellArraysOfSameLength',ZName,GName);
        for k = 1:R
            GNamek = [GName,'{',num2str(k),'}'];
            G{k} = LinearMixedModel.validateGroupingVar(G{k},GNamek,N);
        end
        
    end % validateXYZG.
                    
    function Nsim = validateNsim(Nsim)
%validateNsim - Validates the 'Nsim' input argument to compare method.
%   Nsim = validateNsim(Nsim) takes a potential value of Nsim and returns
%   the validated value. The validated value is always an integer.
%
%   What is checked?
%
%   (1) Either Nsim is empty 
%         or
%   (2) Nsim is a scalar integer >= 0.
%
%   If (1) or (2) does not hold, an error message is thrown. If Nsim is
%   empty, the output Nsim is equal to 0.

        % (1) If Nsim = [], then set Nsim = 0.
        if isempty(Nsim)
            Nsim = 0;
        end
        
        % (2) Ensure that Nsim is a scalar integer >= 0.        
        assertThat(internal.stats.isScalarInt(Nsim,0,Inf),'stats:LinearMixedModel:MustBeNonNegativeInteger','Nsim');

    end % end of validateNsim.
        
    function verbose = validateVerbose(verbose)
%validateVerbose - Validate the parameter 'Verbose'.
%   verbose = validateVerbose(verbose) takes a potential value verbose for
%   the 'Verbose' flag and validates it. If not valid, an error message is
%   thrown.
%
%   What is checked?
%
%   (1) verbose must be a scalar logical (true or false).
%
%   If (1) is not satisfied, then an error message is thrown.
                        
        verbose = ...
            LinearMixedModel.validateLogicalScalar(verbose,'stats:LinearMixedModel:BadVerbose');

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
            LinearMixedModel.validateLogicalScalar(checkhessian,'stats:LinearMixedModel:BadCheckHessian');

    end % end of validateCheckHessian.
    
end

% Public, static fitting methods.
methods(Static, Access='public', Hidden=true)
    
    function model = fit(ds,formula,pmax,varargin)
%   The FIT method is not intended to be called directly. Use FITLME to 
%   create a LinearMixedModel by fitting to data.
%
%   See also FITLME.
        
        if isa(ds,'dataset')
            ds = dataset2table(ds);
        end

        % (1) Ensure that ds and formula are sensible.            
            assertThat(isa(ds,'table'),'stats:LinearMixedModel:Fit_firstinput');
                                 
            assertThat(internal.stats.isString(formula),'stats:LinearMixedModel:Fit_secondinput');
            formula = convertStringsToChars(formula);
            [varargin{:}] = convertStringsToChars(varargin{:});
        
        % (2) Create an empty LME model.
        model = nirs.math.LinearMixedModel();
        
        % (3) Try to parse the formula using all variables in ds.
        model.Formula = classreg.regr.LinearMixedFormula(formula,ds.Properties.VariableNames);                        
        
        % (4) How many grouping variable? How many total observations?
            R = length(model.Formula.RELinearFormula);
            N = size(ds,1);
            
        % (5) Initialize default values for optional parameters.
        
            % 5(a) dfltCovariancePattern
                dfltCovariancePattern = cell(R,1);
                % Default is 'Full' with Cholesky for each grouping variable.
                dfltCovariancePattern(1:R) = ...
                    {classreg.regr.lmeutils.covmats.CovarianceMatrix.TYPE_FULLCHOLESKY}; 

            % 5(b) dfltFitMethod
                % Default is maximum likelihood.
                dfltFitMethod = 'ML'; 

            % 5(c) dfltWeights and dfltExclude
                % Default is ones(N,1), N = size(ds,1).
                dfltWeights = ones(N,1);   
                % Default is to use all (non-NaN) obs.
                dfltExclude = false(N,1);  

            % 5(d) dfltDummyVarCoding
            dfltDummyVarCoding = 'reference';

            % 5(e) dfltOptimizer
            dfltOptimizer = 'quasinewton';
            
            % 5(f) dfltOptimizerOptions. Do not change this. We use an
            % empty value of optimizeroptions to indicate that no options
            % were supplied in the validateOptimizerAndOptions method.
            dfltOptimizerOptions = [];
            
            % 5(g) dfltStartMethod
            dfltStartMethod = 'default';
            
            % 5(h) dfltVerbose
            dfltVerbose = false;
        
            % 5(i) dfltCheckHessian
            dfltCheckHessian = false;
            
        % (6) Process optional parameter name/value pairs.        
        paramNames =   {'CovariancePattern',   'FitMethod',   'Weights',   'Exclude',   'DummyVarCoding',   'Optimizer',   'OptimizerOptions',   'StartMethod',   'Verbose',   'CheckHessian'};
        paramDflts = {dfltCovariancePattern, dfltFitMethod, dfltWeights, dfltExclude, dfltDummyVarCoding, dfltOptimizer, dfltOptimizerOptions, dfltStartMethod, dfltVerbose, dfltCheckHessian};        
        [covariancepattern,fitmethod,weights,exclude,dummyvarcoding,optimizer,optimizeroptions,startmethod,verbose,checkhessian,setflag] = internal.stats.parseArgs(paramNames, paramDflts, varargin{:});                                
        
        % (7) Validate parameter values except covariancepattern.
        fitmethod      = LinearMixedModel.validateFitMethod(fitmethod);
        weights        = LinearMixedModel.validateWeights(weights,N);
        exclude        = LinearMixedModel.validateExclude(exclude,N);
        dummyvarcoding = LinearMixedModel.validateDummyVarCoding(dummyvarcoding);         
        [optimizer,optimizeroptions] ...
                       = LinearMixedModel.validateOptimizerAndOptions(optimizer,optimizeroptions);                     
           startmethod = LinearMixedModel.validateStartMethod(startmethod);
               verbose = LinearMixedModel.validateVerbose(verbose);
          checkhessian = LinearMixedModel.validateCheckHessian(checkhessian);
          
        % (8) verbose will override the Display field in optimizeroptions.
        if (setflag.Verbose == true)
            % User has supplied a 'Verbose' input.
            if (verbose == true)
                optimizeroptions.Display = 'iter';
            else
                optimizeroptions.Display = 'off';
            end
        end
               
        % (9) At this point, optimizer and optimizeroptions should be
        % synchronized with each other and optimizeroptions should be
        % completely filled out with the default options overridden by the
        % user supplied options.        
        
        % (10) Store fitmethod, dummyvarcoding and optimization options in 
        % FitMethod, DummyVarCoding, Optimizer, OptimizerOptions,
        % StartMethod and CheckHessian properties. weights and exclude are 
        % stored in ObservationInfo and covariancepattern is stored in the 
        % object after validation later.
            model.FitMethod = fitmethod;
            model.DummyVarCoding = dummyvarcoding;            
            model.Optimizer = optimizer;
            model.OptimizerOptions = optimizeroptions;
            model.StartMethod = startmethod;
            model.CheckHessian = checkhessian;
            
        % (11) Provisionally fill out VariableInfo and ObservationInfo. We
        % will alter PredLocs and RespLocs later in selectVariables. This 
        % also indirectly validates weights and exclude arguments.
            model.PredictorTypes = 'mixed'; % Grouping vars are predictors.
            model = assignData(model,ds,[],weights,[],...
                model.Formula.VariableNames,exclude);          
        
        % (12) Select our variables and observations. After this point
        % model.ObservationInfo.Subset marks the subset of data to be used
        % for fitting after removing excluded and missing values.
            model = selectVariables(model);
            model = selectObservations(model,exclude);                       
                      
        % (13) Get y, FixedInfo, RandomInfo and GroupingInfo for the 
        % observations that are going to be used in the fit. Error if
        % subset is all false.
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
            
        % (14) Now that we know the size of various random effects design
        % matrices, validate covariancepattern.
        covariancepattern = LinearMixedModel.validateCovariancePattern...
            (covariancepattern,R,model.RandomInfo.q);

        % (15) Store covariancepattern in CovariancePattern property.
        model.CovariancePattern = covariancepattern;                       
        
        % (16) Fit the model by calling StandardLinearMixedModel in fitter.        
        %   doFit is a FitObject method that does the following:
        %   (1) model = selectVariables(model);            % Set PredLocs, RespLoc and update VariableInfo.InModel
        %   (2) model = selectObservations(model,exclude); % Update ObservationInfo.Missing, .Excluded and .Subset
        %   (3) model = fitter(model);                     % Do the actual fitting.
        %   (4) model = postFit(model);                    % Do post fitting.           
        model.pmax=pmax;
        model = doFit(model);

        % (17) Omit excluded points from range  .
        model = updateVarRange(model);               
        
    end % end of fit.
    
    function model = fitmatrix(X,Y,Z,G,varargin)
%   The FITMATRIX method is not intended to be called directly. Use FITLME 
%   or FITLMEMATRIX to create a LinearMixedModel by fitting to data.
%
%   See also FITLME, FITLMEMATRIX.
        [varargin{:}] = convertStringsToChars(varargin{:});
        % (1) Validate matrix inputs X, Y, Z and G.
        [X,Y,Z,G] = LinearMixedModel.validateXYZG(X,Y,Z,G,'X','Y','Z','G');
            
        % (2) Get dimensions.
        p = size(X,2);
        R = length(Z);
        q = zeros(R,1);
        for k = 1:R
            q(k) = size(Z{k},2);
        end
        
        % (3) Default names for X, Y, Z and G.
        % (3a) X.
        dfltFixedEffectPredictors = internal.stats.numberedNames('x',1:p);
        % (3b) Y.
        dfltResponseVarName = 'y';
        % (3c) Z.
        dfltRandomEffectPredictors = cell(R,1);
        for k = 1:R
            zk = ['z',num2str(k)];
            dfltRandomEffectPredictors{k} = ...
                internal.stats.numberedNames(zk,1:q(k));
        end
        % (3d) G.
        dfltRandomEffectGroups = internal.stats.numberedNames('g',1:R);
        
        % (4) Process optional parameter name/value pairs.        
        paramNames =   {'FixedEffectPredictors',   ...
            'RandomEffectPredictors',...
            'ResponseVarName',...
            'RandomEffectGroups'};
        paramDflts = {dfltFixedEffectPredictors,...
            dfltRandomEffectPredictors,...
            dfltResponseVarName,...
            dfltRandomEffectGroups};
        [fepredictors,repredictors,respvar,regroups,~,otherArgs] = ...
            internal.stats.parseArgs(paramNames, paramDflts,varargin{:});
        
        % (5) Validate names for X, Y, Z and G.
        % (5a) Names for X.
        fepredictors = ...
            LinearMixedModel.validateCellVectorOfStrings(fepredictors,...
            'FixedEffectPredictors',p,true);        
        % (5b) Name for Y.
        respvar = ...
            LinearMixedModel.validateString(respvar,'ResponseVarName');        
        % (5c) Names for Z. First ensure that repredictors is a cell
        % array and then ensure that each element of repredictors is a
        % cell vector of strings of the right size.
        repredictors = ...
            LinearMixedModel.validateCellVector(repredictors,...
            'RandomEffectPredictors',R);
        for k = 1:R            
            repredictorskname = ...
                getString(message('stats:LinearMixedModel:String_repredictors_k',num2str(k)));
            repredictors{k} ...
                = LinearMixedModel.validateCellVectorOfStrings(...
                repredictors{k},repredictorskname,q(k),true);
        end        
        % (5d) Names for G. G need not have unique strings.
        regroups = ...
            LinearMixedModel.validateCellVectorOfStrings(regroups,...
            'RandomEffectGroups',R,false);        
        
        % (6) At this point, X, Y, Z, G and their column names have been
        % validated. Now transform X, Y, Z, G and their column names into a
        % table and a formula string.
        [ds,formula] = LinearMixedModel.convertXYZGToDataset(X,Y,Z,G,...
            fepredictors,respvar,repredictors,regroups);
                                
        % (7) Now call the LinearMixedModel.fit method with formula input.
        % Pass in otherArgs{:} to .fit as varargin.
        model = LinearMixedModel.fit(ds,formula,otherArgs{:});
        
        % (8) Store original names for columns of X, Y, Z and G. 
        model.XYZGNames.XNames = fepredictors;
        model.XYZGNames.YName  = respvar;
        model.XYZGNames.ZNames = repredictors;
        model.XYZGNames.GNames = regroups;
        
    end % end of fitmatrix.
    
end % end of methods(Static, Access='public', Hidden).

% Public methods of LinearMixedModel.
methods(Access='public')
    
    function [D,gnames] = designMatrix(model,designtype,gnumbers)
%designMatrix Extracts the fixed or random effects design matrices.
%   X = designMatrix(LME) or designMatrix(LME,'Fixed') returns the N-by-P
%   fixed effects design matrix X for the linear mixed effects model LME
%   where N is the number of observations and P is the number of fixed
%   effects.
%
%   D = designMatrix(LME,'Random') returns the overall random effects
%   design matrix corresponding to a vector B of all random effects in the
%   linear mixed effects model LME. Suppose LME has R grouping variables
%   named g_1,...,g_R. Let Q_1,...,Q_R be the length of random effects
%   vectors associated with g_1,...,g_R respectively. Also, suppose
%   g_1,...,g_R have levels M_1,...,M_R respectively. Then B will be a
%   column vector of length Q_1*M_1 + Q_2*M_2 + ... + Q_R*M_R. B is made by
%   concatenating the BLUPs of random effects vectors corresponding to each
%   level of each grouping variable in the following order:
%
%   g_1 level 1, g_1 level 2,..., g_1 level M_1 followed by
%   g_2 level 1, g_2 level 2,..., g_2 level M_2 followed by
%      .            .                .           
%   g_R level 1, g_R level 2,..., g_R level M_R
%
%   If LME has N observations then D will be of size N-by-length(B) such
%   that D*B is a N-by-1 vector that represents the contribution of all
%   random effects to the response of the LME.
%
%   DSUB = designMatrix(LME,'Random',GNUMBERS) returns a submatrix of the
%   full random effects design matrix. GNUMBERS is a length K integer array
%   with elements in the range 1 to R. DSUB is a subset of the full random
%   effects design matrix corresponding to the grouping variable names
%   indicated by integers in GNUMBERS. For example, suppose GNUMBERS is
%   [1,R] then this specifies only grouping variables g_1 and g_R. Let BSUB
%   be a vector made by concatenating BLUPs of random effects vectors
%   corresponding to each level of g_1 and g_R in the following order:
%
%   g_1 level 1, g_1 level 2,..., g_1 level M_1 followed by         
%   g_R level 1, g_R level 2,..., g_R level M_R
%
%   Then DSUB will be a N-by-length(BSUB) matrix such that DSUB*BSUB
%   represents the contribution of all random effects corresponding to
%   grouping variables g_1 and g_R to the response of the LME. If GNUMBERS
%   is empty, the full random effects design matrix is returned.
%
%   [DSUB,GNAMES] = designMatrix(LME,DESIGNTYPE,GNUMBERS) also returns a 
%   K-by-1 cell array containing the names of grouping variables
%   corresponding to integers in GNUMBERS if DESIGNTYPE is 'Random'. If
%   DESIGNTYPE is 'Fixed' then GNAMES is [] and GNUMBERS is ignored.
%        
%   Example: Fit a quadratic model with continuous and categorical fixed
%            effects. Look at a few rows of the design matrix showing the
%            constant term, the Weight term, and dummy variables for the
%            Origin term.
%      load carsmall
%      ds = dataset(MPG,Weight,Model_Year,Origin);
%      lme = fitlme(ds,'MPG ~ Weight + Origin + (1|Model_Year)');
%      dm = designMatrix(lme);
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
        
    end % end of design matrix.
            
    function [beta,betanames,fetable] = fixedEffects(model,varargin)
%fixedEffects Returns estimates of fixed effects and related statistics.
%   BETA = fixedEffects(LME) returns a vector of estimated fixed effects
%   from a fitted linear mixed effects model LME.
%
%   [BETA,BETANAMES] = fixedEffects(LME) also returns a dataset array
%   BETANAMES containing the name of each fixed effects coefficient in
%   BETA.
%
%   [BETA,BETANAMES,STATS] = fixedEffects(LME) also returns a dataset array
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
%   [BETA,BETANAMES,STATS] = fixedEffects(LME,'PARAM','VALUE',...) also
%   specifies optional parameter name/value pairs to control the fixed
%   effects statistics:
%  
%           'Name'     'Value'
%       'DFMethod'     Specifies the method to use for computing the
%                      approximate degrees of freedom (DF) for the t
%                      statistics that test the fixed effects coefficients
%                      against 0. Options are 'Satterthwaite','Residual'
%                      and 'None'. If 'DFMethod' is 'Satterthwaite', a
%                      Satterthwaite approximation is used to compute DF.
%                      If 'DFMethod' is 'Residual', the DF values are
%                      assumed to be constant and equal to (N-P) where N is
%                      the number of observations and P is the number of
%                      fixed effects. If 'DFMethod' is 'None', then all DF
%                      values are set to infinity. Default is 'Residual'.
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
%      ds = dataset(MPG,Weight,Model_Year);
%      lme = fitlme(ds,'MPG ~ Weight + (1|Model_Year)')
%      fixedEffects(lme)
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
%   B = randomEffects(LME) returns estimates of the best linear unbiased
%   predictors (BLUPs) of all random effects in the linear mixed effects
%   model LME. Suppose LME has R grouping variables named g_1,...,g_R. Let
%   Q_1,...,Q_R be the length of random effects vectors associated with
%   g_1,...,g_R respectively. Also, suppose g_1,...,g_R have levels
%   M_1,...,M_R respectively. Then B will be a column vector of length
%   Q_1*M_1 + Q_2*M_2 + ... + Q_R*M_R. B is made by concatenating the BLUPs
%   of random effects vectors corresponding to each level of each grouping
%   variable in the following order:
%
%   g_1 level 1, g_1 level 2,..., g_1 level M_1 followed by
%   g_2 level 1, g_2 level 2,..., g_2 level M_2 followed by
%      .            .                .           
%   g_R level 1, g_R level 2,..., g_R level M_R
%
%   [B,BNAMES] = randomEffects(LME) also returns a table array BNAMES
%   containing the names of the coefficients in B.
%
%   [B,BNAMES,STATS] = randomEffects(LME) also returns a table array
%   STATS containing the estimates of random effects and related
%   statistics. Table STATS has one row for each random effect associated
%   with a particular grouping variable, level and predictor name with the
%   following columns:
%
%       Group       Grouping variable associated with this random effect
%       Level       Level within the grouping variable
%       Name        Name of the random effect coefficient
%       Estimate    BLUP of random effect
%       SEPred      Standard error of (BLUP minus the random effect)
%       tStat       t statistic for a test that the random effect is zero
%       DF          Estimated degrees of freedom for the t statistic
%       pValue      p-value for the t statistic
%       Lower       Lower limit of a 95% confidence interval
%       Upper       Upper limit of a 95% confidence interval
%
%   [B,BNAMES,STATS] = randomEffects(LME,'PARAM','VALUE',...) specifies
%   optional parameter name/value pairs to control the calculation of
%   random effects statistics:
%  
%           'Name'     'Value'
%       'DFMethod'     Specifies the method to use for computing the
%                      approximate degrees of freedom (DF) for the t
%                      statistics that test the random effects coefficients
%                      against 0. Options are 'Satterthwaite','Residual'
%                      and 'None'. If 'DFMethod' is 'Satterthwaite', a
%                      Satterthwaite approximation is used to compute DF.
%                      If 'DFMethod' is 'Residual', the DF values are
%                      assumed to be constant and equal to (N-P) where N is
%                      the number of observations and P is the number of
%                      fixed effects. If 'DFMethod' is 'None', then all DF
%                      values are set to infinity. Default is 'Residual'.
%
%          'Alpha'     ALPHA, a number between 0 and 1. The columns 'Lower'
%                      and 'Upper' in table STATS will correspond to the
%                      lower and upper limits respectively of 100*(1-ALPHA)
%                      confidence intervals for random effects. Default is
%                      ALPHA=0.05 for 95% confidence intervals.
%
%   Example: Fit a model with a fixed effect and random effects. Get the
%            estimated random effects. They are highly negatively
%            correlated. Verify this by plotting them.
%      load carsmall
%      ds = dataset(MPG,Weight,Model_Year);
%      lme = fitlme(ds,'MPG ~ Weight + (Weight|Model_Year)')
%      re = randomEffects(lme)
%      plot(re(1:2:end),re(2:2:end),'rs')
%
%   See also coefCI, coefTest, fixedEffects.

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
%   PSI = covarianceParameters(LME) extracts the estimated covariance
%   parameters that parameterize the prior covariance of random effects. If
%   the linear mixed effects model LME has R grouping variables named
%   g_1,...,g_R then the output PSI will be a R-by-1 cell array such that
%   PSI{i} will contain the covariance matrix of random effects associated
%   with grouping variable g_i. The order in which grouping variables are
%   assigned numbers 1 to R is the same order in which grouping variables
%   are entered into the FITLME or FITLMEMATRIX functions.
%
%   [PSI,MSE] = covarianceParameters(LME) also extracts an estimate of the
%   residual variance.
%
%   [PSI,MSE,STATS] = covarianceParameters(LME) also returns a cell array
%   STATS of length (R+1) containing covariance parameters and related
%   statistics. STATS{i} is a table array containing statistics on
%   covariance parameters for the i-th grouping variable. STATS{R+1}
%   contains statistics on the residual standard deviation. STATS{i} 
%   contains columns that name each covariance parameter as well as the 
%   following columns:
%
%       Group       Grouping variable name
%       Estimate    Estimate of the covariance parameter
%       Lower       Lower limit of a 95% confidence interval
%       Upper       Upper limit of a 95% confidence interval
%
%   NOTE: It is recommended that the presence or absence of covariance
%   parameters in LME be tested using the COMPARE method, which uses a
%   likelihood ratio test.
%
%   [PSI,MSE,STATS] = covarianceParameters(LME,PARAM1,VALUE1,...)
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
%      ds = dataset(MPG,Weight,Model_Year);
%      lme = fitlme(ds,'MPG ~ Weight + (Weight|Model_Year)')
%      V = covarianceParameters(lme);
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
%FITTED Returns the fitted response from the LME.
%   YFIT = FITTED(LME) returns the N-by-1 vector representing the fitted
%   conditional response from the LME, where N is the number of
%   observations. If the LME has a N-by-P fixed effects design matrix X, an
%   estimated P-by-1 fixed effects vector BETA, a N-by-Q overall random
%   effects design matrix Z and an estimated BLUP B of the overall random
%   effects vector, then YFIT = X*BETA + Z*B.
%
%   YFIT = FITTED(LME,PARAM1,VALUE1,...) accepts optional name/value
%   pairs as follows:
%
%           'Name'     'Value'
%     'Conditional'     Either true or false.  If false then YFIT = X*BETA. 
%                       Default is true.
%
%   Example: Fit a model and plot the residuals vs. fitted values, grouped by
%            Origin.
%      load carsmall
%      ds = dataset(MPG,Weight,Model_Year);
%      lme = fitlme(ds,'MPG ~ Weight + (1|Model_Year)');
%      r = residuals(lme);
%      f = fitted(lme);
%      gscatter(f,r,Origin)
%
%   See also RESIDUALS, RESPONSE, designMatrix.
        [varargin{:}] = convertStringsToChars(varargin{:});
        % (1) Default parameter values.
         dfltConditional = true;

        % (2) Optional parameter names and their default values.
        paramNames =   {'Conditional'};
        paramDflts = {dfltConditional};
           
        % (3) Parse optional parameter name/value pairs.
        wantconditional = ...
            internal.stats.parseArgs(paramNames,paramDflts,varargin{:});   

        % (4) Validate optional parameter values.
        wantconditional ...
            = LinearMixedModel.validateConditional(wantconditional);

        % (5) Forward the call to StandardLinearMixedModel object slme.
        yr = fitted(model.slme,wantconditional);

        % (6) Divide each element of yr by sqrt(w) where w = observation 
        % weights.        
            % (1) First, get the weights for observations used in fit.
            reduce = true;            
            w = getCombinedWeights(model,reduce);
            % (2) Now divide yr(i) by sqrt(w(i)).
            yr = yr./sqrt(w);        
                
        % (7) Fill in NaNs in the output for observations not used in fit.
        % TO DO: LinearModel does not fill in NaNs for missing obs.
        subset = model.ObservationInfo.Subset;
        yfit = NaN(length(subset),1);
        yfit(subset) = yr;  
        
    end % end of fitted.

    function res = residuals(model,varargin)
%RESIDUALS Returns the residuals from a fitted model.
%   R = RESIDUALS(LME) returns the N-by-1 vector of raw conditional
%   residuals from a fitted linear mixed effects model LME where N is the
%   number of observations.
%
%   R = RESIDUALS(LME,PARAM1,VALUE1,...) accepts optional name/value
%   pairs as follows:
%
%           'Name'     'Value'
%     'Conditional'     Either true or false. If false then marginal 
%                       residuals are returned. Default is true. 
%
%    'ResidualType'     Valid values are 'Raw', 'Standardized' and 
%                       'Pearson'. Default is 'Raw'.
%
%   Suppose the LME has a N-by-1 response vector Y, a N-by-P fixed effects
%   design matrix X, an estimated P-by-1 fixed effects vector BETA_HAT, an 
%   N-by-Q overall random effects design matrix Z and an estimated Q-by-1
%   random effects vector B_HAT. Define,
%
%               RC = Y - (X*BETA_HAT + Z*B_HAT)
%           RCTRUE = Y - (X*BETA + Z*B)
%
%               RM = Y - X*BETA_HAT
%           RMTRUE = Y - X*BETA
%
%   Here BETA is the true fixed effects vector and B is the random effects
%   vector. B_HAT is the estimated Best Linear Unbiased Predictor (BLUP) of
%   B. Here's how the elements of residual vector are calculated for each
%   observation i and for various combinations of values of 'Conditional'
%   and 'ResidualType'. STD in the table below, stands for the estimated
%   standard deviation of the quantity in parentheses.
%
%     ResidualType          Conditional=true         Conditional=false
%     ------------          ----------------         -----------------
%        'Raw'                   RC(i)                    RM(i)
%    'Standardized'         RC(i)/STD(RC(i))          RM(i)/STD(RM(i))
%      'Pearson'          RC(i)/STD(RCTRUE(i))      RM(i)/STD(RMTRUE(i))
%
%   Example: Fit a model and plot the residuals vs. fitted values, grouped by
%            Origin.
%      load carsmall
%      ds = dataset(MPG,Weight,Model_Year);
%      lme = fitlme(ds,'MPG ~ Weight + (1|Model_Year)');
%      r = residuals(lme);
%      f = fitted(lme);
%      gscatter(f,r,Origin)
%
%   See also FITTED, RESPONSE, designMatrix.
        [varargin{:}] = convertStringsToChars(varargin{:});
        % (1) Default parameter values.
         dfltConditional = true;
        dfltResidualType = 'Raw';

        % (2) Optional parameter names and their default values.
        paramNames =   {'Conditional',   'ResidualType'};
        paramDflts = {dfltConditional, dfltResidualType};
           
        % (3) Parse optional parameter name/value pairs.
        [wantconditional,residualtype] ...
            = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});   

        % (4) Validate optional parameter values.
        wantconditional ...
            = LinearMixedModel.validateConditional(wantconditional);
        residualtype = LinearMixedModel.validateResidualType(residualtype);
        
        % (5) Forward the call to member variable slme.
        r = residuals(model.slme,wantconditional,residualtype);

        % (6) If we are returning raw residuals, divide each element of r
        % by sqrt(w) where w = observation weights.
        % TODO: LinearModel applies the weight scaling below to all
        % residual types. Need to do the same here.
        if strcmpi(residualtype,'Raw')
            % (1) First, get the weights for observations used in fit.
            reduce = true;            
            w = getCombinedWeights(model,reduce);
            % (2) Now divide r(i) by sqrt(w(i)).
            r = r./sqrt(w);
        end
        
        % (7) Fill in NaNs in the output for observations not used in fit.
        subset = model.ObservationInfo.Subset;
        res = NaN(length(subset),1);
        res(subset) = r;        
        
    end % end of residuals.
    
    function [table,siminfo] = compare(model,altmodel,varargin)
%COMPARE Compare fitted linear mixed effects models.
%   TABLE = COMPARE(LME,ALTLME) performs a likelihood ratio test comparing
%   linear mixed effects models LME and ALTLME which have been fit to the
%   same response vector but with different model specifications. LME must
%   be nested within ALTLME, i.e., LME is obtained from ALTLME by setting
%   some parameters in ALTLME to fixed values (such as 0). The output TABLE
%   is a table array containing the results of the likelihood ratio test. 
%   TABLE has 2 rows, the first row is for LME, and the second row is for 
%   ALTLME. TABLE has the following columns:
%
%           Model       The name of the model.
%           DF          The number of free parameters in the model.
%           AIC         Akaike information criterion for the model.
%           BIC         Bayesian information criterion for the model.
%           LogLik      The maximized log-likelihood for the model.
%           LRStat      Likelihood ratio test statistic for comparing 
%                       ALTLME versus LME.
%           deltaDF     DF for ALTLME minus DF for LME.
%           pValue      p-value for the likelihood ratio test.
%
%   The null and alternative hypotheses are as follows:
%
%       H0: Observed response vector was generated by model LME. 
%       H1: Observed response vector was generated by model ALTLME.
%
%   A p-value for the likelihood ratio test is computed by comparing the
%   observed likelihood ratio test statistic with a chi-squared reference
%   distribution with degrees of freedom deltaDF. A small p-value (e.g., 
%   < 0.05) leads to a rejection of H0 in favor of H1 and acceptance of 
%   model ALTLME. On the other hand, a large p-value (e.g., > 0.05) 
%   reflects insufficient evidence to accept model ALTLME.
%
%       (1) It is recommended that LME and ALTLME be fitted using maximum
%       likelihood (ML) prior to model comparison. If restricted maximum
%       likelihood (REML) is used in fitting LME and ALTLME, then LME must
%       be nested in ALTLME, and in addition LME and ALTLME must have
%       exactly the same fixed effects design matrix.
%  
%       (2) p-values computed using the chi-squared reference distribution
%       as described above can be:
%            - conservative, when testing for the presence or absence of
%              random effects terms and
%            - anti-conservative, when testing for the presence or absence
%              of fixed effects terms
%
%       Hence, testing for fixed effects should be performed either using
%       method fixedEffects or the simulated likelihood ratio test.
%
%   [TABLE,SIMINFO] = COMPARE(LME, ALTLME,'NSim',NSIM,'alpha',ALPHA)
%   performs a simulated likelihood ratio test comparing linear mixed
%   effects models LME and ALTLME which have been fit to the same response
%   vector but with different model specifications.
%
%          'NAME'         'VALUE'
%          'NSim'         NSIM, an integer that specifies the number of 
%                         simulations to use in the simulated likelihood 
%                         ratio test. You must specify the value of 'Nsim'
%                         to do a simulated likelihood ratio test. Default
%                         is 0 for a theoretical likelihood ratio test.
% 
%         'Alpha'         ALPHA, a value between 0 and 1 to specify a
%                         100*(1-ALPHA)% confidence interval for the 
%                         estimated p-value. 'Alpha' is optional with a 
%                         default value of ALPHA = 0.05.
%   
%   The output TABLE has the following format for a simulated likelihood
%   ratio test:
%
%           Model       The name of the model.
%           DF          The number of free parameters in the model.
%           AIC         Akaike information criterion for the model.
%           BIC         Bayesian information criterion for the model.
%           LogLik      The maximized log-likelihood for the model.
%           LRStat      Likelihood ratio test statistic for comparing 
%                       ALTLME versus LME.
%           pValue      Estimated p-value for the simulated likelihood ratio test.
%           Lower       Lower limit of a confidence interval for true pValue.
%           Upper       Upper limit of a confidence interval for true pValue.
%
%   The output SIMINFO is a structure containing simulation output with the 
%   following fields:
%
%          'NAME'        'DESCRIPTION'
%           nsim        Value set for 'NSim'.
%          alpha        Value set for 'Alpha'.
%      pvalueSim        Simulation based p-value.
%    pvalueSimCI        A 1-by-2 vector containing the (1-ALPHA) 
%                       confidence interval for pvalueSim.
%        deltaDF        The number of free parameters in ALTLME minus the
%                       number of free parameters in LME. 
%            TH0        A vector of simulated likelihood ratio test 
%                       statistics under H0.
%
%   The reference distribution of the likelihood ratio test statistic under
%   H0 is generated by simulation as follows:
% 
%         (1) Generate random data YSIM from the fitted model LME. 
%         (2) Fit the model specified in LME and ALTLME to the simulated 
%             data YSIM.
%         (3) Calculate the likelihood ratio test statistic using results
%             from (2).
%         (4) Repeat steps (1)-(3) NSIM times.
%         (5) This produces a reference distribution of the likelihood 
%             ratio test statistic under H0.
%
%   A p-value for the simulated likelihood ratio test is computed by
%   comparing the observed likelihood ratio test statistic (say T) with 
%   the simulated reference distribution in (5) (stored in vector TH0) as
%   follows:
%
%      p-value = ( (number of values in TH0 >= T) + 1 )/(NSIM + 1)
%
%   To account for the uncertainity in the simulated reference
%   distribution, a 100*(1-ALPHA)% confidence interval for the p-value is
%   also computed. When using the simulated likelihood ratio test:
%
%         (1) It is not required that LME be nested within ALTLME.
%
%         (2) LME and ALTLME can be fitted using either maximum likelihood
%         (ML) or restricted maximum likelihood (REML).
%
%   Hence, the simulated likelihood ratio test can be used to compare
%   arbitrary linear mixed effects models.
%
%   [TABLE,SIMINFO] = COMPARE(LME,ALTLME,...,'PARAM',VALUE,...) specifies
%   additional name/value pairs for the likelihood ratio test:
%
%              'NAME'         'VALUE'
%           'Options'         OPTIONS, a structure that contains options
%                             specifying whether to perform computations
%                             for simulated likelihood ratio test in
%                             parallel, and specifying how to use random
%                             numbers during the simulation. OPTIONS uses
%                             the following fields:
%
%                             'UseParallel', 'UseSubstreams' and 'Streams'
%
%                             For detailed information on these fields type
%                             "help parallelstats". Only valid for a 
%                             simulated likelihood ratio test. Default is 
%                             statset('UseParallel',false).
%
%      'CheckNesting'         Either true or false. If 'CheckNesting' is 
%                             true then we attempt to check if the smaller
%                             model LME is nested in the bigger model
%                             ALTLME. It is necessary that LME be nested in
%                             ALTLME for the theoretical likelihood ratio
%                             test to be valid. An error message is thrown
%                             if the nesting requirement is not satisfied.
%                             Valid for both simulated and theoretical
%                             likelihood ratio test, although the nesting
%                             requirements are much weaker for the
%                             simulated likelihood ratio test. If 'REML' is 
%                             used in fitting LME and ALTLME, then LME and 
%                             ALTLME must have identical fixed effects 
%                             design matrices. Default is false.
%
%   Example: Model gas mileage as a function of car weight, with a random
%            effect due to model year. Compare a model having random
%            intercepts with one also having random slopes.
%      load carsmall
%      ds = dataset(MPG,Weight,Model_Year);
%      lme = fitlme(ds,'MPG ~ Weight + (1|Model_Year)')
%      lme2 = fitlme(ds,'MPG ~ Weight + (1 | Model_Year) + (-1 + Weight|Model_Year)')
%      compare(lme,lme2)
%
%   See also fixedEffects, ANOVA, randomEffects, covarianceParameters.
        [varargin{:}] = convertStringsToChars(varargin{:});
        % (1) Ensure that model and altmodel are LinearMixedModel objects.
           model = LinearMixedModel.validateObjectClass(model,...
               'LME','LinearMixedModel');
        altmodel = LinearMixedModel.validateObjectClass(altmodel,...
            'ALTLME','LinearMixedModel');

        % (2) Default parameter values.
                    dfltNsim = 0;
                   dfltAlpha = 0.05;
                 dfltOptions = statset('UseParallel',false);
            dfltCheckNesting = false;

        % (3) Optional parameter names and their default values.
        paramNames = {  'Nsim',   'Alpha',   'Options',   'CheckNesting'};
        paramDflts = {dfltNsim, dfltAlpha, dfltOptions, dfltCheckNesting};
           
        % (4) Parse optional parameter name/value pairs.
        [nsim,alpha,options,checknesting] ...
            = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});   

        % (5) Validate optional parameter values.
                    nsim = LinearMixedModel.validateNsim(nsim);
                   alpha = LinearMixedModel.validateAlpha(alpha);
                 options = LinearMixedModel.validateOptions(options);
            checknesting = ...
                LinearMixedModel.validateCheckNesting(checknesting);
            
        % (6) Now that options is valid, use it along with dfltOptions to
        % create a full user supplied structure.
        options = statset(dfltOptions,options);                 
        
        % (7) Get names for model and altmodel.
            modelName = inputname(1);
            altmodelName = inputname(2);
            if isempty(modelName) || isempty(altmodelName)
                modelName = 'LME';
                altmodelName = 'ALTLME';
            end
            
        % (8) Decide what to do based on the value of nsim.
        if ( nsim == 0 )
            % Standard LRT.                          
                
            % (1) Ensure that model is "nested" in altmodel. This is not
            % foolproof but at least it is some protection.
            if ( checknesting == true )
                isSimulatedTest = false;
                LinearMixedModel.checkNestingRequirement(model,altmodel,...
                    modelName,altmodelName,isSimulatedTest);
            end
                
            % (2) Get the standard LRT table.
            table = LinearMixedModel.standardLRT(model,altmodel,...
                modelName,altmodelName);
            
            % (3) siminfo is not applicable and set to [].
            siminfo = [];                
                
        else
            % Simulated LRT.
            
            % (1) For simulated LRT, nesting requirements are less strict.
            if ( checknesting == true )
                isSimulatedTest = true;
                LinearMixedModel.checkNestingRequirement(model,altmodel,...
                    modelName,altmodelName,isSimulatedTest);
            end
                
            % (2) Get the simulated LRT table.
            [table,siminfo] = LinearMixedModel.simulatedLRT(model,...
                altmodel,modelName,altmodelName,nsim,alpha,options);
            
        end      

    end % end of compare
    
    function hout = plotResiduals(model,plottype,varargin)
% plotResiduals Plot residuals of fitted model
%   plotResiduals(LME,PLOTTYPE) plots the raw conditional residuals from
%   model LME in a plot of type PLOTTYPE. Valid values for PLOTTYPE are:
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
%    'ResidualType'     Valid values are 'Raw' (default), 'Standardized'
%                       and 'Pearson'.
%
%   For more details on the various residual types, please see the help for
%   the residuals method.
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
%      ds = dataset(MPG,Weight,Model_Year,'ObsNames',obsname);
%      lme = fitlme(ds,'MPG ~ Weight + (1|Model_Year)')
%      plotResiduals(lme,'probability')
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
%   STATS = ANOVA(LME) tests the significance of each fixed effect term in
%   the linear mixed effects model LME and returns a table array STATS.
%   Each fixed effect term reported in STATS is either a continuous
%   variable, a grouping variable or an interaction between two or more
%   variables (continuous or grouping). For each fixed effect term, ANOVA
%   performs an F test (marginal test) that all coefficients representing
%   the fixed effect term are zero. There is one row in STATS for each
%   fixed effect term and the following columns:
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
%   STATS = ANOVA(LME,PARAM1,VALUE1,...) accepts additional parameter
%   name/value pairs:
%  
%           'Name'     'Value'
%       'DFMethod'     Specifies the method to use for computing the
%                      approximate denominator degrees of freedom DF2 for
%                      the F-statistics reported in STATS. Options are
%                      'Satterthwaite','Residual' and 'None'. If 'DFMethod'
%                      is 'Satterthwaite', a Satterthwaite approximation is
%                      used to compute DF2. If 'DFMethod' is 'Residual',
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
%      ds = dataset(MPG,Weight,Model_Year,Cylinders);
%      ds.Cylinders = nominal(ds.Cylinders);
%      lme = fitlme(ds,'MPG ~ Weight + Cylinders + (1|Model_Year)')
%      anova(lme)
%
%   See also coefTest, coefCI, fixedEffects.
        
        stats = anova@classreg.regr.LinearLikeMixedModel(model,varargin{:});
        
    end % end of anova.
    
    function Y = response(model)
%RESPONSE Returns the response vector for the linear mixed effects model.
%   Y = RESPONSE(LME) returns the N-by-1 response vector Y used to fit the
%   linear mixed effects model LME, where N is the number of observations.
%
%   See also FITTED, RESIDUALS, designMatrix.

        % (1) Get the subset of observations in fit.
        subset = model.ObservationInfo.Subset;
        
        % (2) Initialize Y to NaN(N,1) where N = length(subset). Then fill
        % in the response values at the observations used in the fit.
        N = length(subset);
        Y = NaN(N,1);
        Y(subset) = model.y;                

    end % end of response.
               
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
            error(message('stats:LinearMixedModel:BadMsgID',msgID));
        end        
        % (2) Create and throw an MException.
        ME = MException( msg.Identifier, getString(msg) );
        throwAsCaller(ME);        
    end

end % end of assertThat.




 function Xw = scaleMatrixUsingWeights(X,w)
%scaleMatrixUsingWeights - Multiplies each row of X using sqrt(w).
%   Xw = scaleMatrixUsingWeights(X,w) takes a N-by-P matrix X and a N-by-1 
%   weight vector w and outputs a scaled matrix Xw such that: 
%           Xw(i,:) = X(i,:) * sqrt(w(i)).       
        
        % (1) Ensure that X is a numeric, real matrix.
        assert( isnumeric(X) & isreal(X) & ismatrix(X) );
        
        % (2) Get the size of X.
        N = size(X,1);
        
        % (3) Ensure that w is a weight vector of length N.
        w = LinearMixedModel.validateWeights(w,N);

        % (4) Form Xw. If all weights are 1 just copy X into Xw.
        if all(w == 1)
            Xw = X;
        else
            Xw = bsxfun(@times,X,sqrt(w));
        end

    end % end of scaleMatrixUsingWeights.
    
    
    
function out = myFilter( f, y )

if(isempty(y))
    out=y;
    return
end

    % here we are just making the first value zero before filtering to
    % avoid weird effects introduced by zero padding
    y1 = y(1,:);
    
    y = bsxfun(@minus,y,y1);
    
    out = filter(f, 1, y);
    out = bsxfun(@plus,out,sum(f)*y1); % add the corrected offset back

end



function w = wfun(r)
s = mad(r, 0) / 0.6745;
r = r/s/4.685;

w = (1 - r.^2) .* (r < 1 & r > -1);
end

