classdef (AllowedSubclasses = {?classreg.regr.ParametricRegression}) Predictor < classreg.regr.FitObject
%Predictor Fitted predictive regression model.
%   Predictor is an abstract class representing a fitted regression model
%   for predicting a response as a function of predictor variables.
%   You cannot create instances of this class directly.  You must create
%   a derived class by calling the fit method of a derived class such as
%   LinearModel, GeneralizedLinearModel, or NonLinearModel.
%
%   See also LinearModel, GeneralizedLinearModel, NonLinearModel.

%   Copyright 2011-2014 The MathWorks, Inc.

    properties(Dependent,GetAccess='public',SetAccess='protected',Abstract=true)
%Fitted - Vector of fitted (predicted) values.
%   The Fitted property is a vector providing the response values predicted
%   by the regression model. The fitted values are computed using the
%   predictor values used to fit the model. Use the PREDICT method to
%   compute predictions for other predictor values.
%
%   See also Predictor, Residuals, predict, random.
        Fitted

%Residuals - Residual values.
%   The Residuals property is a table of residuals. It is a table array that
%   has one row for each observation and, depending on the type of
%   regression model, has one or more of the following columns: 
%      Raw           Observed minus predicted response values
%      Pearson       Raw residuals normalized by RMSE
%      Standardized  Internally Studentized residuals based on RMSE
%      Studentized   Externally Studentized residuals based on S2_i
%
%   To obtain any of these columns as a vector, index into the property
%   using dot notation. For example, in the LinearModel LM the ordinary or
%   raw residual vector is
%      r = LM.Residuals.Raw
%
%   See also Predictor, Fitted, predict, random.
        Residuals
    end
    
    methods(Abstract, Access='public')
        ypred = predict(model,varargin)
        ysim = random(model,varargin)
    end
    
    methods(Abstract, Access='protected')
        % The predict method would normally take a dataset/table array or a matrix
        % containing all variables.  This method exists to allow prediction
        % with a matrix that contains only the required predictor variables
        % without blowing it up to contain all the variables only to then pull
        % it apart to get the design matrix.
        ypred = predictPredictorMatrix(model,Xpred);
    end
    
    methods(Hidden, Access='public')
        function yPred = predictGrid(model,varargin)
            % Prediction over a grid, works only for numeric predictors
            gridMatrices = gridVectors2gridMatrices(model,varargin);
            outSize = size(gridMatrices{1});
            gridCols = cellfun(@(x) x(:), gridMatrices,'UniformOutput',false);
            predData = table(gridCols{:},'VariableNames',model.PredictorNames);
            yPred = predict(model,predData);
            yPred = reshape(yPred,outSize);
        end
    end
       
    methods(Access='public')
        function yPred = feval(model,varargin)
%FEVAL Evaluate model as a function
%    YPRED = FEVAL(M,X1,X2,...Xp) computes an array YPRED of predicted
%    values from the regression model M using predictor values X1, X2, ...,
%    Xp, where P is the number of predictors in M.  Each X argument must be
%    the same size, and must have the same type as the corresponding
%    predictor variable. The size of the YPRED output is the same as the
%    common size of the X inputs.
%
%    YPRED = FEVAL(M,DS) or YPRED = FEVAL(M,X) accepts a dataset/table DS or
%    matrix X containing values of all of the predictors.
%
%    The PREDICT method can compute confidence bounds as well as predicted
%    values. The M.Fitted property provides predicted values using the
%    predictor values used to fit M.
%
%    Example:
%       % Fit model to car data; superimpose fitted cuves on scatter plot
%       load carsmall
%       d = dataset(MPG,Weight);
%       d.Year = ordinal(Model_Year);
%       lm = fitlm(d,'MPG ~ Year + Weight + Weight^2')
%       w = linspace(min(d.Weight),max(d.Weight))';
%       gscatter(d.Weight, d.MPG, d.Year);
%       line(w, feval(lm,w,'70'), 'Color','r')
%       line(w, feval(lm,w,'76'), 'Color','g')
%       line(w, feval(lm,w,'82'), 'Color','b')
%       
%    See also Fitted, predict, LinearModel, GeneralizedLinearModel, NonlinearModel.

            npreds = model.NumPredictors;
            if isa(varargin{1},'dataset')
                varargin{1} = dataset2table(varargin{1});
            end
            if nargin-1 == npreds && ...% separate predictor variable arguments
               ~(nargin==2 && isa(varargin{1},'table'))
                predArgs = varargin;
                
                % Get common arg length considering possible scalar
                % expansion
                sizeOut = [1 1];
                for i = 1:length(predArgs)
                    thisarg = predArgs{i};
                    if ischar(thisarg)
                        if size(thisarg,1)~=1
                            sizeOut = [size(thisarg,1),1];
                            break
                        end
                    else
                        if ~isscalar(thisarg)
                            sizeOut = size(thisarg);
                            break
                        end
                    end
                end
                
                % Get args as cols, expanding scalars as needed
                asCols = predArgs;
                for i = 1:length(predArgs)
                    thisarg = predArgs{i};
                    if ischar(thisarg)
                        thisarg = cellstr(thisarg);
                    end
                    if isscalar(thisarg)
                        thisarg = repmat(thisarg,sizeOut);
                    elseif ~isequal(size(predArgs{i}),sizeOut)
                        error(message('stats:classreg:regr:Predictor:InputSizeMismatch'));
                    end
                    asCols{i} = thisarg(:);
                end
                
                % Evaluate model on predictors as columns, then resize
                Xpred = table(asCols{:},'VariableNames',model.PredictorNames);
                yPred = reshape(predict(model,Xpred),sizeOut);
            elseif nargin == 2 
                predVars = varargin{1};
                if isa(predVars,'table')
                    yPred = predict(model,predVars);
                else % single predictor variable matrix
                    if size(predVars,2) ~= npreds
                        error(message('stats:classreg:regr:Predictor:BadNumColumns', npreds));
                    end
                    yPred = predictPredictorMatrix(model,predVars);
                end
            else
                error(message('stats:classreg:regr:Predictor:BadNumInputs', npreds, npreds));
            end
        end
    end
    
    methods(Abstract,Access='protected')
        D = get_diagnostics(model,type)
    end
    methods(Access='protected')
        function r = get_residuals(model)
            r = getResponse(model) - predict(model);
        end
        
        function yfit = get_fitted(model)
            yfit = predict(model);
        end
    end
end
