function model = stepwise(X,varargin) % [X, y | DS], start, ...
% Not intended to be called directly. Use STEPWISELM to fit a LinearModel
% using stepwise regression.
%
%   See also STEPWISELM.

%   Copyright 2011-2014 The MathWorks, Inc.

[X,y,haveDataset,otherArgs] = LinearModel.handleDataArgs(X,varargin{:});


paramNames = {'Intercept' 'PredictorVars' 'ResponseVar' 'Weights' 'Exclude' 'CategoricalVars' ...
    'VarNames' 'Lower' 'Upper' 'Criterion' 'PEnter' 'PRemove' 'NSteps' 'Verbose'};
paramDflts = {true [] [] [] [] [] [] 'constant' 'interactions' 'SSE' [] [] Inf 1};

% Default model is constant only.
if isempty(otherArgs)
    start = 'constant';
else
    arg1 = otherArgs{1};
    if mod(length(otherArgs),2)==1 % odd, model followed by pairs
        start = arg1;
        otherArgs(1) = [];
    elseif internal.stats.isString(arg1) && ...
            any(strncmpi(arg1,paramNames,length(arg1)))
        % omitted model but included name/value pairs
        start = 'constant';
    end
end

[intercept,predictorVars,responseVar,weights,exclude,asCatVar, ...
    varNames,lower,upper,crit,penter,premove,nsteps,verbose,supplied] = ...
    internal.stats.parseArgs(paramNames, paramDflts, otherArgs{:});

[penter,premove] = classreg.regr.TermsRegression.getDefaultThresholds(crit,penter,premove);

if ~isscalar(verbose) || ~ismember(verbose,0:2)
    error(message('stats:LinearModel:BadVerbose'));
end

% *** need more reconciliation among start, lower, upper, and between those
% and intercept, predictorVars, varNames

if ~supplied.ResponseVar && (classreg.regr.LinearFormula.isTermsMatrix(start) || classreg.regr.LinearFormula.isModelAlias(start))
    if isa(lower,'classreg.regr.LinearFormula')
        responseVar = lower.ResponseName;
        supplied.ResponseVar = true;
    else
        if internal.stats.isString(lower) && ~classreg.regr.LinearFormula.isModelAlias(lower)
            lower = LinearModel.createFormula(supplied,lower,X, ...
                        predictorVars,responseVar,intercept,varNames,haveDataset);
            responseVar = lower.ResponseName;
            supplied.ResponseVar = true;
        elseif isa(upper,'classreg.regr.LinearFormula')
            responseVar = upper.ResponseName;
            supplied.ResponseVar = true;
        else
            if internal.stats.isString(upper) && ~classreg.regr.LinearFormula.isModelAlias(upper)
                upper = LinearModel.createFormula(supplied,upper,X, ...
                            predictorVars,responseVar,intercept,varNames,haveDataset);
                responseVar = upper.ResponseName;
                supplied.ResponseVar = true;
            end
        end
    end
end
    
if ~isa(start,'classreg.regr.LinearFormula')
    ismodelalias = classreg.regr.LinearFormula.isModelAlias(start);
    start = LinearModel.createFormula(supplied,start,X, ...
                predictorVars,responseVar,intercept,varNames,haveDataset);
else
    ismodelalias = false;
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

% Remove categorical powers, if any, but warning only if these powers were
% requested explicitly, not just created via something like 'quadratic'
nvars = size(X,2);
if haveDataset
    isCat = varfun(@internal.stats.isDiscreteVar,X,'OutputFormat','uniform');
else
    isCat = [repmat(internal.stats.isDiscreteVar(X),1,nvars) internal.stats.isDiscreteVar(y)];
    nvars = nvars+1;
end
if ~isempty(asCatVar)
    isCat = classreg.regr.FitObject.checkAsCat(isCat,asCatVar,nvars,haveDataset,start.VariableNames);
end
if any(isCat)
    start = removeCategoricalPowers(start,isCat,ismodelalias);
    lower = removeCategoricalPowers(lower,isCat,ismodelalias);
    upper = removeCategoricalPowers(upper,isCat,ismodelalias);
end

if haveDataset
    model = LinearModel.fit(X,start.Terms,'ResponseVar',start.ResponseName, ...
        'Weights',weights,'Exclude',exclude,'CategoricalVars',asCatVar,'RankWarn',false);
else
    model = LinearModel.fit(X,y,start.Terms,'ResponseVar',start.ResponseName, ...
        'Weights',weights,'Exclude',exclude,'CategoricalVars',asCatVar, ...
        'VarNames',start.VariableNames,'RankWarn',false);
end

model.Steps.Start = start;
model.Steps.Lower = lower;
model.Steps.Upper = upper;
model.Steps.Criterion = crit;
model.Steps.PEnter = penter;
model.Steps.PRemove = premove;
model.Steps.History = [];

model = stepwiseFitter(model,nsteps,verbose);
checkDesignRank(model);
end
