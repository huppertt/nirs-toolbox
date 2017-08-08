classdef (AllowedSubclasses = {?classreg.regr.Predictor}) FitObject < classreg.learning.internal.DisallowVectorOps & internal.matlab.variableeditor.VariableEditorPropertyProvider
%FitObject Fitted statistical regression.
%   FitObject is an abstract class representing a fitted regression model.
%   You cannot create instances of this class directly.  You must create
%   a derived class by calling the fit method of a derived class such as
%   LinearModel, GeneralizedLinearModel, or NonLinearModel.
%
%   See also LinearModel, GeneralizedLinearModel, NonLinearModel.

%   Copyright 2011-2014 The MathWorks, Inc.


    properties(GetAccess='public', SetAccess='protected')
%VariableInfo - Information about variables used in the fit.
%   The VariableInfo property is a table providing information about the
%   variables used in the fit. There is one row for each variable. The
%   table contains the following columns:
%
%     Class - the class of the variable ('double', 'cell', 'nominal',
%             etc.).
%     Range - the range of values taken by the variable, either a
%             two-element vector of the form [min,max] for a numeric
%             variable, or a cell or categorical array listing all unique
%             values of the variable for a cell or categorical variable.
%     InModel - the logical value true if the variable is a predictor in
%             the fitted model, or false if it is not.
%     IsCategorical - the logical value true if the variable has a type
%             that is treated as a categorical predictor (such as cell,
%             logical, or categorical) or if it is specified as categorical
%             by the 'Categorical' option to the fit method, or false if is
%             treated as a continuous predictor.
%
%   See also FitObject.
        VariableInfo = table({'double'},{[NaN; NaN]},false,false, ...
                               'VariableNames',{'Class' 'Range' 'InModel'  'IsCategorical'}, ...
                               'RowNames',{'y'});

%ObservationInfo - Information about observations used in the fit.
%   The ObservationInfo property is a table providing information about
%   the observations used in the fit. There is one row for each
%   observation. The table contains the following columns:
%
%     Weights - the value of the weight variable for that observation
%             (default 1).
%     Excluded - the logical value true if the observation was excluded
%             from the fit using the 'Exclude' argument, or false
%             otherwise.
%     Missing - the logical value true if the observation was excluded from
%             the fit because any response or predictor value is missing.
%             Missing values include NaN for numeric variables, empty
%             cells for cell arrays, blank rows for char arrays, and the
%             <undefined> value for categorical arrays.
%     Subset - the logical value true if the observation was used in the
%             fit, or false if it was not used because it is missing or
%             excluded.
%
%   See also FitObject.
        ObservationInfo = table(zeros(0,1),false(0,1),false(0,1),false(0,1), ...
                                  'VariableNames',{'Weights' 'Excluded' 'Missing' 'Subset'});
    end
    properties(GetAccess='protected', SetAccess='protected')
        PredLocs = zeros(1,0);
        RespLoc = 1;
        Data = []; % dataset/table(y,x1,x2,...) or struct('X',[],'y',[])
        IsFitFromData = false;
        PredictorTypes = 'numeric'; % or 'mixed'
        ResponseType = 'numeric'; % or 'categorical'
        
        WorkingValues = struct; % scalar struct
    end
    properties(Dependent, GetAccess='public', SetAccess='protected')
%Variables - Table of variables used in the fit.
%   The Variables property is a table containing data for the variables
%   used in the fit. This is a copy of the original dataset/table when the fit is
%   based on a dataset/table, or a dataset/table created from the predictor matrix or 
%   matrices and response vector when the fit is based on those arrays.
%
%   See also FitObject.
        Variables

%NumVariables - Number of variables used in the fit.
%   The NumVariables property gives the number of variables used in the
%   fit. This is the number of variables in the original dataset/table when the
%   fit is based on a dataset/table, or the total number of columns in the
%   predictor matrix or matrices and response vector when the fit is based
%   on those arrays. It includes variables, if any, that are not used as
%   predictors or as the response.
%
%   See also FitObject.
        NumVariables

%VariableNames - Names of variables used in the fit.
%   The VariableNames property is a cell array of strings containing the
%   names of the variables in the fit.
%
%   If the fit is based on a dataset/table, this property provides the names of
%   the variables in that dataset/table.
%
%   If the fit is based on a predictor matrix or matrices and response
%   vector, this property is taken from the variable names supplied to the
%   fitting function (if any). Otherwise the variables are given default 
%   names.
%
%   See also FitObject.
        VariableNames

%NumPredictors - Number of predictors used in the fit.
%   The NumPredictors property gives the number of variables used as
%   predictors in the fit.
%
%   See also FitObject.
        NumPredictors

%PredictorNames - Names of predictors used in the fit.
%   The PredictorNames property is a cell array of strings containing the
%   names of the variables used as predictors in the fit.
%
%   See also FitObject.
        PredictorNames

%ResponseName - Names of the response.
%   The ResponseName property is a character string giving the name of the
%   variable used as the response in the fit.
%
%   See also FitObject.
        ResponseName

%NumObservations - Number of observations used in the fit.
%   The NumObservations property gives the number of observations used in
%   the fit. This is the number of observations supplied in the original
%   dataset/table or matrix, minus any excluded rows or rows with missing (NaN)
%   values.
%
%   See also FitObject.
        NumObservations

%ObservationNames - Names of observations used in the fit.
%   The ObservationNames property is a cell array of strings containing the
%   names of the observations used in the fit. If the fit is based on a
%   dataset/table containing observation names, this property is taken from those
%   names. If the fit is based on a matrix or on a dataset/table without
%   observation names, this property is an empty cell array.
%
%   See also FitObject.
        ObservationNames
    end
    
    methods % get/set methods
        function vars = get.Variables(model)
            vars = getVariables(model); % to allow override
        end
        function n = get.NumVariables(model)
            n = size(model.VariableInfo,1);
        end
        function vnames = get.VariableNames(model)
            vnames = model.VariableInfo.Properties.RowNames;
        end
        function n = get.NumPredictors(model)
            n = length(model.PredictorNames);
        end
        function vnames = get.PredictorNames(model)
            vnames = model.VariableInfo.Properties.RowNames(model.PredLocs);
        end
        function vname = get.ResponseName(model)
            vname = model.VariableInfo.Properties.RowNames{model.RespLoc};
        end
        function n = get.NumObservations(model)
            n = sum(model.ObservationInfo.Subset);
        end
        function onames = get.ObservationNames(model)
            if fitFromDataset(model)
                onames = model.Data.Properties.RowNames; % possibly empty names
            else
                onames = {};
            end
        end
    end % get/set methods
    
    methods(Abstract, Access='protected')
        model = fitter(model);
    end
    methods(Abstract, Hidden, Access='public')
        t = title(model);
    end
    methods(Abstract, Static, Access='public')
        model = fit(varargin);
    end
    methods(Abstract, Hidden, Access='public')
        disp(model)
        val = feval(model,varargin); % requires a modified fnchk if hidden, g701463
    end
    
    methods(Access='protected')
        function model = noFit(model,varNames)
            p = length(varNames) - 1;
            model.PredLocs = 1:p;
            model.RespLoc = p+1;
            viName = varNames;
            viClass = repmat({'double'},p+1,1);
            viRange = repmat({[NaN; NaN]},p+1,1);
            viInModel = [true(p,1); false];
            viIsCategorical = false(p+1,1);
            model.VariableInfo = table(viClass,viRange,viInModel,viIsCategorical, ...
                                     'VariableNames',{'Class' 'Range' 'InModel' 'IsCategorical'}, ...
                                     'RowNames',viName);
            model.Data = table([{zeros(0,p+1)},varNames]);
            model.IsFitFromData = false;
        end
        
        % --------------------------------------------------------------------
        function model = doFit(model,exclude)
            if nargin < 2
                exclude = model.ObservationInfo.Excluded;
            end
            
            % Choose variables before removing missing values
            model = selectVariables(model);
            model = selectObservations(model,exclude);
            
            model = fitter(model);
            model = postFit(model);
            model.IsFitFromData = true;
            
            % Clear out the WorkingValues structure, anything saved in there
            % by subclasses during fitting will have to be recreated when needed
            model.WorkingValues = [];
        end
                
        % --------------------------------------------------------------------
        function model = assignData(model,X,y,w,asCat,varNames,excl)
            % This is a generic implementation suitable for column-oriented
            % numeric/categorical data.  The data can be
            %    * in one dataset/table array X containing matrix or vector variables (and y is []),
            %    * in separate matrices X and y
            % We don't check types, we leave that for the selectVariables method.
            % varNames is required for data in X,y, ignored for data in a dataset/table array
            if isa(X,'dataset')
                X = dataset2table(X);
            end
            haveDataset = isa(X,'table');

            % Do some preliminary processing of the input data
            if haveDataset && all(varfun(@ismatrix,X,'OutputFormat','uniform'))
                [nobs,nvars] = size(X);
                predLocs = 1:(nvars-1);
                respLoc = nvars;
            elseif ismatrix(X) && ismatrix(y)
                % recognize vectors of either orientation, this assumes n > 1
                if isvector(X)
                    X = X(:);  % force to a column
                end
                [nobs,p] = size(X); % p is number of predictors
                if ischar(X)
                    p = 1;
                end
                if isvector(y)
                    y = y(:);  % force to a column
                end
                if size(y,1) ~= nobs
                    error(message('stats:classreg:regr:FitObject:PredictorResponseMismatch'));
                end
                nvars = p + 1;
                predLocs = 1:(nvars-1);
                respLoc = nvars;
            else
                error(message('stats:classreg:regr:FitObject:MatricesRequired'));
            end

            % Process the exclusion vector
            [excl,predLocs] = getDataVariable(excl,nobs,X,predLocs,respLoc,'exclude');
            if islogical(excl)
                if numel(excl)~=nobs
                    error(message('stats:classreg:regr:FitObject:BadExcludeLength',nobs));
                end
            else
                if any(excl<0) || any(excl>nobs) || any(excl~=round(excl))
                    error(message('stats:classreg:regr:FitObject:BadExcludeValues'));
                end
                tmp = excl;
                excl = false(nobs,1);
                excl(tmp) = true;
            end
            if all(excl)
                error(message('stats:classreg:regr:FitObject:AllExcluded'));
            end
            
            % Now compute variable ranges and other information
            if haveDataset
                viName = X.Properties.VariableNames;
                viClass = varfun(@class,X,'OutputFormat','cell');
                
                % Provisionally choose the last variable as the response, and
                % the rest as predictors.
                viInModel = [true(1,nvars-1) false];
                
                viIsCategorical = varfun(@internal.stats.isDiscreteVar,X,'OutputFormat','uniform');

                if ~isempty(asCat)
                    viIsCategorical = classreg.regr.FitObject.checkAsCat(viIsCategorical,asCat,nvars,true,viName);
                end
                viRange = cell(nvars,1);
                for i = 1:nvars
                    viRange{i} = getVarRange(X.(viName{i}),viIsCategorical(i),excl);
                end

                data = X;
                obsNames = X.Properties.RowNames;   
            elseif ismatrix(X) && ismatrix(y)
                % The Variables property puts the response last, match that
                viName = varNames;
                viClass = [repmat({class(X)},1,p) {class(y)}];
                
                % Provisionally choose all variables in X as predictors.
                viInModel = [true(1,p) false];
                
                viIsCategorical = [repmat(internal.stats.isDiscreteVar(X),1,p) internal.stats.isDiscreteVar(y)];
                if ~isempty(asCat)
                    viIsCategorical = classreg.regr.FitObject.checkAsCat(viIsCategorical,asCat,nvars,false,viName);
                end
                viRange = cell(nvars,1);
                for i = 1:(nvars-1)
                    viRange{i} = getVarRange(X(:,i),viIsCategorical(i),excl);
                end
                viRange{end} = getVarRange(y,viIsCategorical(end),excl);
                
                data = struct('X',{X},'y',{y});
                obsNames = {};
            end

            [w,predLocs] = getDataVariable(w,nobs,X,predLocs,respLoc,'weights');
            if isempty(w)
                w = ones(nobs,1);
            elseif any(w<0) || numel(w)~=nobs
                error(message('stats:classreg:regr:FitObject:BadWeightValues', nobs));
            end
                
            % We do not check for empty data, a class may want to allow models
            % that are not fit to data, or may have some edge case behavior.
            model.Data = data;
            model.PredLocs = predLocs; % provisionally
            model.RespLoc = respLoc;   % provisionally
            model.VariableInfo = table(viClass(:),viRange(:),viInModel(:),viIsCategorical(:), ...
                                     'VariableNames',{'Class' 'Range' 'InModel' 'IsCategorical'}, ...
                                     'RowNames',viName);
            model.ObservationInfo = table(w,excl,false(nobs,1),false(nobs,1), ...
                                     'VariableNames',{'Weights' 'Excluded' 'Missing' 'Subset'},...
                                     'RowNames',obsNames);
        end

        % ------------------------------------
        function model = updateVarRange(model)
            % Update variable range of model with VariableInfo already
            % filled in. We'll just update the continuous variables, as we
            % need to preserve the original set of categories and their
            % order for categorical variables used as predictors.
            vrange = model.VariableInfo.Range;
            vcat = model.VariableInfo.IsCategorical;
            excl = ~model.ObservationInfo.Subset;
            vclass = ismember(model.VariableInfo.Class, {'nominal','ordinal','categorical'});
            
            if fitFromDataset(model)
                vnames = model.Data.Properties.VariableNames;
                for i = 1:length(vnames)
                    if ~vcat(i)&&~vclass(i)
                        vrange{i} = getVarRange(model.Data.(vnames{i}),vcat(i),excl);
                    end
                end
            else
                nx = size(model.Data.X,2);
                for i = 1:nx
                    if ~vcat(i)
                        vrange{i} = getVarRange(model.Data.X(:,i),vcat(i),excl);
                    end
                end
                if ~vclass(nx+1)
                    vrange{nx+1} = getVarRange(model.Data.y,vcat(nx+1),excl);
                end
            end
            model.VariableInfo.Range = vrange;
        end
        
        % --------------------------------------------------------------------
        function model = selectVariables(model)
            % Subclasses can override to modify predLocs and respLoc before
            % the type/size checks that follow
            model.VariableInfo.InModel(:) = false;
            model.VariableInfo.InModel(model.PredLocs) = true;
            
            haveDataset = fitFromDataset(model);
            if haveDataset
                data = model.Data;
                isNumVar = varfun(@isnumeric,data,'OutputFormat','uniform');
                isNumVec = isNumVar & varfun(@isvector,data,'OutputFormat','uniform');
                isCatVec = varfun(@internal.stats.isDiscreteVec,data,'OutputFormat','uniform');
                switch model.PredictorTypes
                case 'numeric'
                    if ~all(isNumVar(model.PredLocs)) % allow multiple columns
                        error(message('stats:classreg:regr:FitObject:PredictorMatricesRequired'));
                    end
                case 'mixed'
                    if ~all(isNumVec(model.PredLocs) | isCatVec(model.PredLocs))
                        error(message('stats:classreg:regr:FitObject:PredictorMatricesRequired'));
                    end
                otherwise
                    % subclass does its own check
                end
                isNumVecY = isNumVec(model.RespLoc);
                isCatVecY = isCatVec(model.RespLoc);
                
            else % data in separate X and y
                X = model.Data.X;
                isNumVarX = isnumeric(X) || islogical(X);
                isCatVecX = isa(X,'categorical') && isvector(X);
                switch model.PredictorTypes
                case 'numeric'
                    if ~isNumVarX
                        error(message('stats:classreg:regr:FitObject:PredictorMatricesRequired'));
                    end
                case 'mixed'
                    if ~(isNumVarX || isCatVecX)
                        error(message('stats:classreg:regr:FitObject:PredictorMatricesRequired'));
                    end
                otherwise
                    % subclass needs to do its own check
                end
                y = model.Data.y;
                isNumVecY = isnumeric(y) && isvector(y);
                isCatVecY = internal.stats.isDiscreteVec(y);
            end
            
            % Allow logical responses to be treated as numeric
            if isCatVecY && strcmp(model.ResponseType,'numeric') ...
                         && strcmp(model.VariableInfo.Class(model.RespLoc),'logical')
                model.VariableInfo.IsCategorical(model.RespLoc) = false;
                isCatVecY = false;
                isNumVecY = true;
            end
            
            switch model.ResponseType
            case 'numeric'
                if ~isNumVecY
                    error(message('stats:classreg:regr:FitObject:NonNumericResponse'));
                end
            case 'categorical'
                if ~isCatVecY
                    error(message('stats:classreg:regr:FitObject:NonCategoricalResponse'));
                end
            otherwise
                % subclass needs to do its own check
            end
            if model.VariableInfo.IsCategorical(model.RespLoc) && ~strcmp(model.ResponseType,'categorical')
                error(message('stats:classreg:regr:FitObject:ResponseTypeMismatch'));
            end
        end
        
        % --------------------------------------------------------------------
        function model = selectObservations(model,exclude,missing)
            % The third input allows subclasses to override and compute their
            % own logical vector before calling this method.
            
            % Number of observations in the data, not necessarily the number of
            % observations used for the fit.
            nobs = size(model.ObservationInfo,1);
            haveDataset = fitFromDataset(model);
            
            % Default handling for missing values: remove them
            if nargin < 3 || isempty(missing)
                if haveDataset
                    vn = model.VariableNames;
                    y = model.Data.(vn{model.RespLoc});
                    missing = internal.stats.hasMissingVal(y(:,1));
                    for j = model.PredLocs
                        x = model.Data.(vn{j});
                        missing = missing | internal.stats.hasMissingVal(x);
                    end
                elseif isnumeric(model.Data.X)
                    missing = any(isnan(model.Data.X),2) | isnan(model.Data.y(:,1));
                elseif isa(model.Data.X,'categorical')
                    missing = any(isundefined(model.Data.X),2) | isnan(model.Data.y(:,1));
                else
                    missing = isnan(model.Data.y(:,1));
                end
            else
                if ~isvector(missing) || length(missing) ~= nobs
                    error(message('stats:classreg:regr:FitObject:BadMissingLength'));
                end
                missing = missing(:); % force to a column
            end
                        
            if isempty(exclude)
                exclude = false(nobs,1);
            else
                [isObsIndices,isInt] = internal.stats.isIntegerVals(exclude,1,nobs);
                if isObsIndices && isvector(exclude) % observation indices
                    where = exclude;
                    exclude = false(nobs,1);
                    exclude(where) = true;
                elseif isInt && all(exclude>0) && isvector(exclude) % give another msg for binary double vector
                    error(message('stats:classreg:regr:FitObject:BadExcludeIndices'));
                elseif haveDataset && internal.stats.isStrings(exclude) % observation names
                    [tf,where] = ismember(exclude,X.Properties.ObsNames);
                    if ~all(tf)
                        error(message('stats:classreg:regr:FitObject:BadExcludeNames'));
                    end
                    exclude = false(nobs,1);
                    exclude(where) = true;
                elseif islogical(exclude) && isvector(exclude)
                    if length(exclude) ~= nobs
                        error(message('stats:classreg:regr:FitObject:BadExcludeLength',nobs));
                    end
                else
                    error(message('stats:classreg:regr:FitObject:BadExcludeType'));
                end
            end
            exclude = exclude(:); % force to a column
            
            subset = ~(missing | exclude);
            model.ObservationInfo.Missing = missing;
            model.ObservationInfo.Excluded = exclude;
            model.ObservationInfo.Subset = subset;
        end
        
        % --------------------------------------------------------------------
        function model = postFit(model)
        end
        
        % --------------------------------------------------------------------
        function tf = hasData(model)
            tf = ~isempty(model.Data);
        end
        function tf = fitFromDataset(model)
            tf = isa(model.Data,'dataset')||isa(model.Data,'table');
        end     

        % --------------------------------------------------------------------
        function [varName,varNum] = identifyVar(model,var)
            [tf,varName] = internal.stats.isString(var,true); % accept a scalar cellstr
            if tf
                varNum = find(strcmp(varName,model.VariableNames));
                if isempty(varNum)
                    error(message('stats:classreg:regr:FitObject:UnrecognizedName', varName));
                end
            elseif internal.stats.isScalarInt(var,1)
                varNum = var;
                if varNum > model.NumVariables
                    error(message('stats:classreg:regr:FitObject:BadVariableNumber', model.NumVariables));
                end
                varName = model.VariableNames{varNum};
            else
                error(message('stats:classreg:regr:FitObject:BadVariableSpecification'));
            end
        end
        
        % --------------------------------------------------------------------
        function vars = getVariables(model,iobs,varargin) % dependent property get/set method override
            % get a dataset/table containing all or some of the fit's variables,
            % possibly only a subset of the observations.  iobs and jvars can
            % be anything that dataset/table subscripting supports.
            if nargin < 2
                iobs = ':';
                varargin{1} = 1:model.NumVariables;
            end
            if fitFromDataset(model)
                vars = model.Data(iobs,varargin{:});
            else % Data is a struct with fields X and y
                % By convention, put response variable last
                Xvars = num2cell(model.Data.X(iobs,:),1);
                vars = table(Xvars{:},model.Data.y(iobs,:),'VariableNames',model.VariableNames);
                vars = vars(:,varargin{:});
            end
        end
        
        % --------------------------------------------------------------------
        function [varData,varName,varNum] = getVar(model,var)
            % get a single one of the fit's variables.  var can be
            % a var name or an index.
            [varName,varNum] = identifyVar(model,var);
            if fitFromDataset(model)
                varData = model.Data.(varName);
            else % Data is a struct with fields X and y
                if varNum < model.NumVariables
                    varData = model.Data.X(:,varNum);
                else
                    varData = model.Data.y;
                end
            end
            t = ~model.ObservationInfo.Subset;
            if any(t)
                if isnumeric(varData)
                    varData(t) = NaN;
                elseif isa(varData,'categorical')
                    varData(t) = {''};
                end
            end
        end
        
        % --------------------------------------------------------------------
        function var = getResponse(model)
            % get the fit's response variable.
            if fitFromDataset(model)
                var = model.Data.(model.VariableNames{model.RespLoc});             
            else % Data is a struct with fields X and y
                var = model.Data.y;
            end
        end
        % --------------------------------------------------------------------
        function X = getData(model)
            if fitFromDataset(model)
                X = model.Data;
            else % Data is a struct with fields X and y
                X = model.Data.X;
            end
        end
        
        % --------------------------------------------------------------------
        function [Xeval,respSz] = preEval(model,toMatrix,varargin)
            % This is a generic implementation suitable for column-oriented
            % numeric data.  The predictor data can either be in a single
            % numeric matrix or dataset/table array containing all the variables as
            % columns, in separate numeric arrays (one per predictor) with
            % common shape.
            p = model.NumPredictors;
            if isa(varargin{1},'dataset')
                varargin{1} = dataset2table(varargin{1});
            end
            if isscalar(varargin) && isa(varargin{1},'table') % data in a dataset/table array
                Xeval = varargin{1};
                [tf,predLocs] = ismember(model.PredictorNames,Xeval.Properties.VariableNames);
                if ~all(tf)
                    error(message('stats:classreg:regr:FitObject:BadPredictorName'));
                elseif ~reconcilePredictorTypes(model,Xeval)
                    error(message('stats:classreg:regr:FitObject:BadPredictorType'));
                end
                respSz = [size(Xeval,1) 1];
                if toMatrix
                    Xeval = varfun(@(x) x, Xeval,'InputVariables', predLocs);
                    Xeval = cat(2,Xeval{:});
                end
            elseif (nargin-2) == p
                % data in separate arrays
                Xeval = varargin;
                if ~reconcilePredictorTypes(model,Xeval)
                    error(message('stats:classreg:regr:FitObject:BadPredictorType'));
                end
                if p > 1
                    varSz = size(Xeval{1});
                    if ~all(cellfun(@(v) isequal(size(v),varSz), Xeval))
                        error(message('stats:classreg:regr:FitObject:PredictorSizeMismatch'));
                    end
                end
                respSz = varSz;
                if toMatrix
                    Xeval = cat(length(varSz)+1,Xeval{:});
                    Xeval = reshape(Xeval,prod(varSz),p);
                end
            elseif isscalar(varargin)
                Xeval = varargin{1};
                if isnumeric(Xeval) && ismatrix(Xeval) % data in one matrix
                    if size(Xeval,2) ~= p
                        error(message('stats:classreg:regr:FitObject:BadPredictorColumns'));
                    elseif ~reconcilePredictorTypes(model,Xeval)
                        error(message('stats:classreg:regr:FitObject:BadPredictorType'));
                    end
                    respSz = [size(Xeval,1) 1];
                    if ~toMatrix
                        Xeval = num2cell(Xeval,1);
                    end
                else
                    error(message('stats:classreg:regr:FitObject:BadPredictorInput'));
                end
            else
                error(message('stats:classreg:regr:FitObject:BadPredictorCount'))
            end
        end
        
        % --------------------------------------------------------------------
        function tf = reconcilePredictorTypes(model,vars) %#ok<INUSD>
            if hasData(model)
                tf = true;
            else
                tf = true;
            end
        end
                
        % --------------------------------------------------------------------
        function gm = gridVectors2gridMatrices(model,gv)
            p = model.NumPredictors;
            
            if ~(iscell(gv) && isvector(gv) && length(gv) == p)
                error(message('stats:classreg:regr:FitObject:BadGridSize', p));
            elseif ~reconcilePredictorTypes(model,gv)
                error(message('stats:classreg:regr:FitObject:BadGridTypes', p));
            end
            
            if ~all(cellfun(@(x)isnumeric(x)&&isvector(x),gv))
                error(message('stats:classreg:regr:FitObject:NonNumericGrid', p));
            end
            
            if p > 1
                [gm{1:p}] = ndgrid(gv{:});
            else
                gm = gv; % ndgrid would treat this as ndgrid(gv,gv)
            end
        end
    end % protected
        
    methods(Hidden, Access='public')
        % --------------------------------------------------------------------
        function val = fevalgrid(model,varargin)
            gridMatrices = gridVectors2gridMatrices(model,varargin);
            val = feval(model,gridMatrices{:});
        end
        
        % --------------------------------------------------------------------
        function [varargout] = subsref(a,s)
            switch s(1).type
            case '()'
                error(message('stats:classreg:regr:FitObject:ParenthesesNotAllowed'));
            case '{}'
                error(message('stats:classreg:regr:FitObject:NoCellIndexing'));
            case '.'
                % Look ahead so that references such as fit.Variables.x1 do
                % not require creating all of fit.Variables.  Let the
                % built-in handle other properties.
                if strcmp(s(1).subs,'Variables') && ~isscalar(s) && ...
                        ~(isequal(s(2).type,'.') && isequal(s(2).subs,'Properties'))
                    if isequal(s(2).type,'.')
                        p = getVar(a,s(2).subs);
                    else
                        p = getVariables(a,s(2).subs{:});
                    end
                    if length(s) > 2
                        [varargout{1:nargout}] = builtin('subsref',p,s(3:end));
                    else
                        [varargout{1:min(nargout,1)}] = p;
                    end
                else
                    % the builtin subsref handles regular property and
                    % method access, allows access to hidden properties,
                    % and forbids access to protected properties
                    [varargout{1:nargout}] = builtin('subsref',a,s);
                end
            end
        end
        function a = subsasgn(a,s,~)
            switch s(1).type
            case '()'
                error(message('stats:classreg:regr:FitObject:NoParetnthesesAssignment'));
            case '{}'
                error(message('stats:classreg:regr:FitObject:NoCellAssignment'));
            case '.'
                if any(strcmp(s(1).subs,properties(a)))
                    error(message('stats:classreg:regr:FitObject:ReadOnly', s( 1 ).subs, class( a )));
                else
                    error(message('stats:classreg:regr:FitObject:NoMethodProperty', s( 1 ).subs, class( a )));
                end
            end
        end
        
    end % hidden public
    
    methods(Hidden, Static, Access='public')
        function a = empty(varargin) %#ok<STOUT>
            error(message('stats:classreg:regr:FitObject:NoEmptyAllowed'));
        end
        
        function model = loadobj(obj)
            if isa(obj.VariableInfo,'dataset')
                obj.VariableInfo = dataset2table(obj.VariableInfo);
            end
            if isa(obj.ObservationInfo,'dataset')
                obj.ObservationInfo = dataset2table(obj.ObservationInfo);
            end
            if isa(obj.Data,'dataset')
                obj.Data.Properties.Description = '';
                obj.Data = dataset2table(obj.Data);
            end
            model = obj;
        end
        
    end % hidden static public
    
    methods(Static, Access='protected')
        function opts = checkRobust(robustOpts)
            
            if isequal(robustOpts,'off') || isempty(robustOpts) || isequal(robustOpts,false)
                opts = [];
                return;
            end
            if internal.stats.isString(robustOpts) || isa(robustOpts,'function_handle') || isequal(robustOpts,true)
                if isequal(robustOpts,'on') || isequal(robustOpts,true)
                    wfun = 'bisquare';
                else
                    wfun = robustOpts;
                end
                robustOpts = struct('RobustWgtFun',wfun,'Tune',[]);
            end
            if isstruct(robustOpts)
                % A robust structure must have a weight function. It may
                % have a tuning constant. For pre-defined weight functions,
                % will determine the default weighting function later. For
                % function handles, set the value to 1.
                fn = fieldnames(robustOpts);
                if ~ismember('RobustWgtFun',fn) || isempty(robustOpts.RobustWgtFun)
                    opts = [];
                    return
                end
                if ~ismember('Tune',fn)
                    robustOpts.Tune = [];
                end
                if internal.stats.isString(robustOpts.RobustWgtFun)
                    [opts.RobustWgtFun,opts.Tune] = dfswitchyard('statrobustwfun',robustOpts.RobustWgtFun,robustOpts.Tune);
                    if isempty(opts.Tune)
                        error(message('stats:classreg:regr:FitObject:BadRobustName'));
                    end
                else
                    opts = struct('RobustWgtFun',robustOpts.RobustWgtFun,'Tune',robustOpts.Tune);
                end
            else
                error(message('stats:classreg:regr:FitObject:BadRobustValue'));
            end
        end
        % ----------------------------------------------------------------------------
        function [varNames,predictorVars,responseVar] = ...
                getVarNames(varNames,predictorVars,responseVar,nx)
            if isempty(varNames)
                % Create varNames from predictorVars and responseVar, supplied or
                % default
                if ~isempty(predictorVars) && (iscell(predictorVars) || ischar(predictorVars))
                    if iscell(predictorVars)
                        predictorVars = predictorVars(:);
                    else
                        predictorVars = cellstr(predictorVars);
                    end
                    if length(predictorVars)~=nx || ...
                            ~internal.stats.isStrings(predictorVars,true)
                        error(message('stats:classreg:regr:FitObject:BadPredNames'));
                    end
                    pnames = predictorVars;
                else
                    pnames = internal.stats.numberedNames('x',1:nx)';
                    if ~isempty(responseVar) && internal.stats.isString(responseVar)
                        pnames = genvarname(pnames,responseVar);
                    end
                    if isempty(predictorVars)
                        predictorVars = pnames;
                    else
                        predictorVars = pnames(predictorVars);
                    end
                end
                if isempty(responseVar)
                    responseVar = genvarname('y',predictorVars);
                end
                varNames = [pnames; {responseVar}];
            else
                if ~internal.stats.isStrings(varNames,true)
                    error(message('stats:classreg:regr:FitObject:BadVarNames'));
                end
                
                % If varNames is given, figure out the others or make sure they are
                % consistent with one another
                if isempty(responseVar) && isempty(predictorVars)
                    % Default is to use the last name for the response
                    responseVar = varNames{end};
                    predictorVars = varNames(1:end-1);
                    return
                end
                
                % Response var must be a name
                if ~isempty(responseVar)
                    [tf,rname] = internal.stats.isString(responseVar,true);
                    if tf
                        responseVar = rname;
                        if ~ismember(responseVar,varNames)
                            error(message('stats:classreg:regr:FitObject:MissingResponse'));
                        end
                    else
                        error(message('stats:classreg:regr:FitObject:BadResponseVar'))
                    end
                end
                
                % Predictor vars must be names or an index
                if ~isempty(predictorVars)
                    [tf,pcell] = internal.stats.isStrings(predictorVars);
                    if tf
                        predictorVars = pcell;
                        if ~all(ismember(predictorVars,varNames))
                            error(message('stats:classreg:regr:FitObject:InconsistentNames'));
                        end
                    elseif isValidIndexVector(varNames,predictorVars)
                        predictorVars = varNames(predictorVars);
                    else
                        error(message('stats:classreg:regr:FitObject:InconsistentNames'))
                    end
                end
                
                % One may still be empty
                if isempty(predictorVars)
                    predictorVars = setdiff(varNames,{responseVar});
                elseif isempty(responseVar)
                    % If predictorVar is given, there should be just the response left
                    responseVar = setdiff(varNames,predictorVars);
                    if isscalar(responseVar)
                        responseVar = responseVar{1};
                    else
                        error(message('stats:classreg:regr:FitObject:AmbiguousResponse'));
                    end
                else
                    if ~ismember({responseVar},varNames) || ...
                            ~all(ismember(predictorVars,varNames)) || ...
                            ismember({responseVar},predictorVars)
                        error(message('stats:classreg:regr:FitObject:InconsistentNames'))
                    end
                end
            end
        end
        function asCat = checkAsCat(isCat,asCat,nvars,haveDataset,VarNames)
            [isVarIndices,isInt] = internal.stats.isIntegerVals(asCat,1,nvars);
            if isVarIndices && isvector(asCat) % var indices
                where = asCat;
                asCat = false(nvars,1);
                asCat(where) = true;
            elseif isInt && isvector(asCat)
                error(message('stats:classreg:regr:FitObject:BadAsCategoricalIndices'));
            elseif internal.stats.isStrings(asCat) % variable names
                [tf,where] = ismember(asCat,VarNames);
                if ~all(tf)
                    error(message('stats:classreg:regr:FitObject:BadAsCategoricalNames'));
                end
                asCat = false(nvars,1);
                asCat(where) = true;
            elseif islogical(asCat) && isvector(asCat)
                if (haveDataset && length(asCat)~=nvars) || ...
                        (~haveDataset && length(asCat)~=nvars-1)
                    error(message('stats:classreg:regr:FitObject:BadAsCategoricalLength'));
                end
                if ~haveDataset
                    asCat = [asCat(:)',false]; % include response as well
                end
            else
                error(message('stats:classreg:regr:FitObject:BadAsCategoricalType'));
            end
            asCat = asCat(:)';
            asCat = (isCat | asCat);
        end
    end % static protected
end


function range = getVarRange(v,asCat,excl)
v(excl,:) = [];
if asCat % create a list of the unique values
    if isa(v,'categorical')
        % For categorical classes, get the values actually present in the
        % data, not the set of possible values.
        range = unique(v(:));
        range = range(~isundefined(range)); % empty if NaNs
            
    else
        % For classes other than categorical, the list of unique values is
        % also the list of possible values.  But get it in the same order as
        % what grp2idx defines for each class.
        [~,~,range] = grp2idx(v); % leaves NaN or '' out of glevels
            
    end
    if ~ischar(range)
        range = range(:)'; % force a row
    end
elseif isnumeric(v) || islogical(v) % find min and max
    range = [min(v,[],1)  max(v,[],1)]; % ignores NaNs unless all NaNs
else
    range = NaN(1,2);
end
end

function [w,predLocs] = getDataVariable(w,~,X,predLocs,respLoc,vtype)
if isempty(w)
    return
end
if isa(X,'dataset') 
    X = dataset2table(X);
end
if isa(X,'table') && internal.stats.isString(w) % a dataset variable name
    [tf,wloc] = ismember(w,X.Properties.VariableNames);
    if ~tf
        error(message('stats:classreg:regr:FitObject:BadVariableName', vtype, w));
    end
    w = X.(w);
    predLocs = setdiff(predLocs,wloc);
    if wloc == respLoc
        respLoc = max(predLocs);
        predLocs = setdiff(predLocs,respLoc);
    end
end
if ~(isnumeric(w) || islogical(w)) || ~isvector(w) || ~isreal(w)
    error(message('stats:classreg:regr:FitObject:BadVariableValues', vtype));
end
w = w(:); % force to a column
end


function tf = isValidIndexVector(A,idx)
if isempty(idx)
    tf = true;
elseif ~isvector(idx)
    tf = false;
elseif islogical(idx)
    tf = (length(idx) == length(A));
elseif isnumeric(idx)
    tf = all(ismember(idx,1:length(A)));
else
    tf = false;
end
end
