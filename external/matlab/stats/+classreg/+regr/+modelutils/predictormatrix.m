function [X,cols2vars,gnames,glevels] = predictormatrix(data,varargin)
%PREDICTORMATRIX Construct predictor variable matrix from dataset/table array.
%   X = PREDICTORMATRIX(A) returns a numeric predictor variable matrix X
%   created from the variables in the NOBS-by-NVARS dataset/table array A.
%   PREDICTORMATRIX copies numeric variables in A directly to columns of X,
%   and uses GRP2IDX to convert variables in A that are categorical, logical,
%   cell arrays of strings, or char, into group indices and copies those to
%   columns of X.  Only numeric variables in A may have more than one column,
%   and all variables in A must have two dimensions.
%
%   X is NOBS-by-NCOLS, where NCOLS is the sum of the numbers of columns in
%   the predictor variables in A.  X is a double matrix unless any of the
%   predictor variables in A are single.
%
%   PREDICTORMATRIX treats the first through next-to-last variables in A as
%   predictor variables, and creates columns in X corresponding to those
%   variables. PREDICTORMATRIX treates the last variable in A as a response
%   variable, and does not create any corresponding columns in X.
%
%   [X,COLS2VARS] = PREDICTORMATRIX(A) returns an NVARS-by-NCOLS matrix
%   VARS2COLS, whose (I,J)th element contains a 1 if the J-th column of X
%   corresponds to the I-th variable in A, and 0 otherwise.
%
%   [X,COLS2VARS,GNAMES] = PREDICTORMATRIX(A) returns the cell array GNAMES
%   containing the group names of the columns of X that correspond to grouping
%   variables in A.  GNAMES{I} is empty if X(:,I) corresponds to a column in A
%   that is not a grouping variable.
%
%   [X,COLS2VARS,GNAMES,GLEVELS] = PREDICTORMATRIX(A) returns the cell array
%   GLEVELS containing the group levels of the columns of X that correspond to
%   grouping variables in A.  GLEVELS{I} is empty if X(:,I) corresponds to a
%   column in A that is not a grouping variable.
%
%
%   [...] = PREDICTORMATRIX(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%   optional parameter name/value pairs to control how PREDICTORMATRIX creates
%   columns in X from variables in A.  Parameters are:
%
%      'PredictorVars'    Specifies which variables in A to treat as predictors.
%                         Only these variables are copied/converted to columns
%                         in X. Default is 1:NVARS-1.
%      'ResponseVar'      Specifies which variable in A to treat as the
%                         response.  This variable is not copied/converted to
%                         columns in X.  Default is 1.
%      'CategoricalVars'  Specifies which variables in A to treat as categorical
%                         variables.  These variables are converted to group
%                         indices before copying to columns in X.  Default is
%                         false for numeric variables, true for categorical,
%                         logical, cell arrays of strings or char.
%
%   See also DATASET, DESIGNMATRIX, DATASET/DOUBLE, DATASET/SINGLE.

%   Copyright 2010-2014 The MathWorks, Inc.

if isa(data,'dataset')
    data = dataset2table(data);
end

[nobs,nvars] = size(data);
varnames = data.Properties.VariableNames;

catVars = varfun(@(x) isa(x,'categorical') || iscellstr(x) || ischar(x) || islogical(x), data,'OutputFormat','uniform');

paramNames = {'predictorvars' 'responsevar' 'categoricalvars'};
paramDflts = {   1:(nvars-1)         nvars           catVars };
[predictorVars,responseVar,treatAsCategorical,supplied] = ...
    internal.stats.parseArgs(paramNames, paramDflts, varargin{:});

haveDataset = isa(data,'table');

% Default is to treat everything except the response as predictors, unless
% there's a terms matrix.
if supplied.predictorvars
    if isnumeric(predictorVars)
        if any(predictorVars < 1) || any(predictorVars > nvars)
            error(message('stats:classreg:regr:modelutils:BadPredVarsNumeric'));
        end
    elseif islogical(predictorVars)
        if length(predictorVars) ~= nvars
            error(message('stats:classreg:regr:modelutils:BadPredVarsLogical'));
        end
        predictorVars = find(predictorVars);
    elseif haveDataset && internal.stats.isStrings(predictorVars)
        [~,predVarInds] = ismember(predictorVars,varnames);
        if any(predVarInds == 0)
            ind0 = find(predVarInds == 0,1,'first');
            error(message('stats:classreg:regr:modelutils:BadPredVarsName',predictorVars{ind0}));
        end
        predictorVars = predVarInds;
    else
        if haveDataset
            error(message('stats:classreg:regr:modelutils:BadPredVarsDataset'));
        else
            error(message('stats:classreg:regr:modelutils:BadPredVarsMatrix'));
        end
    end
    
    if ~supplied.responsevar
        % Treat the one variable not in the predictors as the response.
        responseVar = setdiff(1:nvars,predictorVars);
        if numel(responseVar) > 1 % none is OK
            error(message('stats:classreg:regr:modelutils:CannotDetermineResponse'));
        end
    end
end

% Default is to treat the first var as the response, unless there's a terms
% matrix.
if supplied.responsevar
    if isempty(responseVar)
        % OK
    elseif isscalar(responseVar) && isnumeric(responseVar)
        if (responseVar < 1) || (responseVar > nvars)
            error(message('stats:classreg:regr:modelutils:BadResponseNumeric'));
        end
    elseif haveDataset && internal.stats.isString(responseVar)
        responseVarIdx = find(strcmp(responseVar,varnames));
        if isempty(responseVarIdx)
            error(message('stats:classreg:regr:modelutils:BadResponseName',responseVar));
        end
        responseVar = responseVarIdx;
    else
        if haveDataset
            error(message('stats:classreg:regr:modelutils:BadResponseDataset'));
        else
            error(message('stats:classreg:regr:modelutils:BadResponseMatrix'));
        end
    end
    
    if ~supplied.predictorvars
        % Treat everything except the response as predictors.
        predictorVars = setdiff(1:nvars,responseVar);
    elseif ~isempty(intersect(responseVar,predictorVars))
        error(message('stats:classreg:regr:modelutils:SamePredictorResponse'));
    end
end

if any(varfun(@ndims,data,'InputVariables',predictorVars,'OutputFormat','uniform') ~= 2)
    error(message('stats:classreg:regr:modelutils:TwoDimVariable'));
end

varClasses = varfun(@class,data,'InputVariables',predictorVars,'OutputFormat','cell');
charVars = strcmp('char',varClasses);

ncols = varfun(@(x)size(x,2),data,'InputVariables',predictorVars,'OutputFormat','uniform');
ncols(charVars) = 1;

if strcmp('single',varClasses)
    outClass = 'single';
else
    outClass = 'double';
end
sumncols = sum(ncols);
X = zeros(nobs,sumncols,outClass);
cols2vars = zeros(nvars,sumncols);
gnames = cell(1,sumncols);
glevels = cell(1,sumncols);
k = 0;
for j = 1:length(predictorVars)
    Xj = data.(varnames{predictorVars(j)});
    kk = k + (1:ncols(j));
    if treatAsCategorical(j)
        [X(:,kk),gnames{kk},glevels{kk}] = grp2idx(Xj); % errors if Xj not a colun vector
    elseif isnumeric(Xj) || islogical(Xj)
        X(:,kk) = Xj;
    else
        if catVars(j)
            error(message('stats:classreg:regr:modelutils:CatNotContinuous'));
        else
            error(message('stats:classreg:regr:modelutils:BadVariableType'));
        end
    end
    cols2vars(j,kk) = 1;
    k = k + ncols(j);
end
