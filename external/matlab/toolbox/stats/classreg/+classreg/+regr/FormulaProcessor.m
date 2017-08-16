classdef FormulaProcessor

%   Copyright 2011-2012 The MathWorks, Inc.

  properties(Abstract,Constant,GetAccess='protected')
        rules
        p
        irules
  end
    properties(Abstract,Access='protected')
        isMultivariate
    end
    properties(Constant,GetAccess='protected')
        allElseStr = '...';
    end
    
    properties(GetAccess='protected',SetAccess='protected')
        tree = [];
        str = 'y ~ 0';
        
        % private copies of dependent public properties
        responseName = 'y';
        varNames = {'y'};
        link = 'identity';
        terms = zeros(0,1);
    end
    
    properties(GetAccess='public',SetAccess='protected')
        FunctionCalls = cell(1,0);
    end

    properties(Dependent,GetAccess='public',SetAccess='public')
        ResponseName
        VariableNames % names for columns of Terms
        Link
        Terms % a terms-by-variables matrix
    end
    properties(Dependent,GetAccess='public',SetAccess='protected')
        ModelFun
        InModel         % which vars are in model as a predictor?
        LinearPredictor % the thing on the RHS of the '~'
        PredictorNames  % names of predictors actually in model
        TermNames       % names of the rows of Terms
        HasIntercept
        NTerms      % number of rows in Terms
        NVars       % number of columns in Terms
        NPredictors % number of predictors actually in model
    end
    methods
        function terms = get.Terms(f)
            terms = f.terms;
        end
        function f = set.Terms(f,update)
            if classreg.regr.FormulaProcessor.isTermsMatrix(update)
                if size(update,2) ~= size(f.terms,2)
                    error(message('stats:classreg:regr:LinearFormula:TermsBadColumns', size( f.terms, 2 )));
                end
                f.terms = removeDupRows(update);
            elseif internal.stats.isString(update)
                [ustr,utree] = parseStr(f,update);
                if utree(1,1) ~= f.irules.LinearPredictor
                    error(message('stats:classreg:regr:LinearFormula:TermsBadLinearPredictor'));
                end
                f.terms = sortTerms(createTerms(f,utree,ustr,1,f.varNames));
            else
                error(message('stats:classreg:regr:LinearFormula:TermsMatrixOrString'));
            end
            f.str = getStr(f);
        end
        function link = get.Link(f)
            link = f.link;
        end
        function f = set.Link(f,update)
            dfswitchyard('stattestlink',update,'double');
            f.link = update;
            f.str = getStr(f);
        end
        function name = get.ResponseName(f)
            name = f.responseName;
        end
        function f = set.ResponseName(f,update)
            if ~internal.stats.isString(update)
                error(message('stats:classreg:regr:LinearFormula:ResponseNameNotString'));
            end
            f.responseName = update;
            f.str = getStr(f);
        end
        function names = get.VariableNames(f)
            names = f.varNames;
        end
        function f = set.VariableNames(f,update)
            [tf,update] = internal.stats.isStrings(update);
            if ~tf || length(update) ~= size(f.terms,2)
                error(message('stats:classreg:regr:LinearFormula:BadVariableNames', size( f.terms, 2 )));
            end
            f.varNames = update;
            f.str = getStr(f);
        end
        
        function fun = get.ModelFun(~)
            fun = @(b,X) X*b;
        end
        function tf = get.InModel(f)
            tf = any(f.terms>0,1);
        end
        function names = get.PredictorNames(f)
            names = f.varNames(f.InModel);
        end
        function names = get.TermNames(f)
            names = classreg.regr.modelutils.terms2names(f.terms,f.varNames); % in term order, not sorted
        end
        function tf = get.HasIntercept(f)
            % When terms has rows and columns, look for an all-zero row.  When
            % it no rows, then there's no intercept.  When it has a row
            % but no columns, it does have an intercept.  (It won't have
            % multiple rows and no columns.)
            tf = any(all(f.terms==0,2),1);
        end
        function expr = get.LinearPredictor(f)
            expr = terms2expr(f.terms,f.varNames);
        end
        function n = get.NVars(f)
            n = size(f.terms,2); % how many vars in the universe
        end
        function n = get.NPredictors(f)
            n = sum(f.InModel); % how many vars appear in Terms
        end
        function n = get.NTerms(f)
            n = size(f.terms,1);
        end
    end
        
    methods(Access='public')
        function f = FormulaProcessor(modelSpec,varNames,responseVar,hasIntercept,link)
            if nargin < 1, return, end
            
            if nargin < 2, varNames = []; end
            if nargin < 3, responseVar = []; end
            if nargin < 4, hasIntercept = []; end
            if nargin < 5, link = []; end
            
            if classreg.regr.FormulaProcessor.isTermsMatrix(modelSpec)
                % The columns in a terms matrix define the universe of
                % variables (possible predictors and response), and these may
                % optionally be named with varNames.  A column of zeros
                % indicates a variable that is not a predictor.  The response
                % variable must have (or is found as) a zero column.
                nvars = size(modelSpec,2);
                if isMissingArg(varNames)
                    % wait until we know the response name to create default var names
                elseif internal.stats.isStrings(varNames,true)
                    if numel(varNames) ~= nvars
                        error(message('stats:classreg:regr:LinearFormula:BadVarNameLength'));
                    end
                    varNames = mustBeUnique(asCellStr(varNames));
                else
                    error(message('stats:classreg:regr:LinearFormula:BadVarNameValue'));
                end
                if isMissingArg(responseVar)
                    respLoc = find(all(modelSpec==0,1));
                    if ~isscalar(respLoc)
                        error(message('stats:classreg:regr:LinearFormula:UndeterminedResponse'));
                    end
                    if isMissingArg(varNames)
                        responseName = 'y'; %#ok<*PROP>
                    else
                        responseName = varNames{respLoc};
                    end
                elseif internal.stats.isString(responseVar)
                    responseName = responseVar;
                    if isMissingArg(varNames)
                        respLoc = find(all(modelSpec==0,1));
                        if ~isscalar(respLoc)
                            error(message('stats:classreg:regr:LinearFormula:UndeterminedResponse'));
                        end
                    else
                        respLoc = find(strcmp(responseName,varNames));
                        if isempty(respLoc)
                            error(message('stats:classreg:regr:LinearFormula:BadResponseName'));
                        elseif any(modelSpec(:,respLoc)~=0)
                            error(message('stats:classreg:regr:LinearFormula:ResponseInTerms'));
                        end
                    end
                elseif internal.stats.isScalarInt(responseVar,1,nvars)
                    respLoc = responseVar;
                    if any(modelSpec(:,respLoc)~=0)
                        error(message('stats:classreg:regr:LinearFormula:ResponseIsPredictor'));
                    end
                    if isMissingArg(varNames)
                        responseName = 'y';
                    else
                        responseName = varNames{respLoc};
                    end
                else
                    error(message('stats:classreg:regr:LinearFormula:ResponseNameOrNumber'));
                end
                if isMissingArg(varNames)
                    varNames = {};
                    varNames([1:(respLoc-1) (respLoc+1):nvars]) = internal.stats.numberedNames('x',1:nvars-1);
                    varNames{respLoc} = responseName;
                end
                if isMissingArg(link)
                    link = 'identity';
                end
                f.terms = removeDupRows(modelSpec);
                f.varNames = varNames;
                f.responseName = responseName;
                if size(f.terms,1) < size(modelSpec,1)
                    warning(message('stats:classreg:regr:LinearFormula:RemoveDupTerms'));
                end
                f.link = link;
                f.FunctionCalls = cell(1,0);
                f.str = composeFormulaString(f.link,f.responseName,terms2expr(f.terms,f.varNames));
                
                % Intercept ordinarily would not be specified, because it's
                % implicit in the terms matrix.  Allow it if no conflict.
                if ~isMissingArg(hasIntercept) && ~isequal(hasIntercept,f.HasIntercept)
                    error(message('stats:classreg:regr:LinearFormula:InterceptFlagConflict'));
                end
            elseif classreg.regr.FormulaProcessor.isModelAlias(modelSpec) ...
                    || (iscell(modelSpec) && classreg.regr.FormulaProcessor.isModelAlias(modelSpec{1}))
                % A model alias requires that you define the universe of vars
                % with varNames, or at least the total number of variables.
                % The responseName may be specified, or is taken to be the
                % last variable.  Ordinarily, all variables except the
                % response will be used as predictors.  You may also specify
                % which variables you actually want as predictors, using
                % {alias,predNames}, in which case you must specify varNames.
                % The response may be specified, or is the one variable not
                % specified as a predictor.
                if isMissingArg(varNames)
                    error(message('stats:classreg:regr:LinearFormula:MissingVarNameNumber'));
                elseif internal.stats.isStrings(varNames)
                    varNames = mustBeUnique(asCellStr(varNames));
                    haveVarNames = true;
                    nvars = numel(varNames);
                elseif internal.stats.isScalarInt(varNames) % number of vars given
                    haveVarNames = false;
                    nvars = varNames;
                    % wait to create default var names until we know the
                    % predictor names or the response name
                else
                    error(message('stats:classreg:regr:LinearFormula:BadVarNames'));
                end
                if isMissingArg(responseVar)
                    % wait to create default response name until we know
                    % whether to get it from the var names or not
                    responseName = [];
                    respLoc = [];
                elseif internal.stats.isString(responseVar)
                    responseName = responseVar;
                    respLoc = [];
                elseif internal.stats.isScalarInt(responseVar,1,nvars)
                    responseName = [];
                    respLoc = responseVar;
                elseif islogical(responseVar) && isvector(responseVar) && ...
                       length(responseVar)<=nvars && ...
                       (sum(responseVar)==1 || f.isMultivariate)
                    responseName = [];
                    respLoc = find(responseVar);
                else
                    error(message('stats:classreg:regr:LinearFormula:InvalidResponseName'));
                end
                if iscell(modelSpec) % modelSpec == {alias,predVars}
                    alias = modelSpec{1};
                    predVars = modelSpec{2};
                    
                    [isVarIndices,isInt] = internal.stats.isIntegerVals(predVars,1,nvars);
                    if isVarIndices && isvector(predVars) % var indices
                        where = predVars;
                        predVars = false(nvars,1); predVars(where) = true;
                    elseif isInt && isvector(predVars)
                        error(message('stats:classreg:regr:LinearFormula:PredictorsOutOfRange'));
                    elseif internal.stats.isStrings(predVars) % variable names
                        if haveVarNames
                            [tf,predLocs] = ismember(predVars,varNames);
                            if ~all(tf)
                                error(message('stats:classreg:regr:LinearFormula:PredictorNotVar'));
                            end
                            predVars = false(1,nvars); predVars(predLocs) = true;
                        else
                            predVars = [true(1,nvars-1) false];
                        end
                    elseif islogical(predVars) && isvector(predVars)
                        if length(predVars) ~= nvars
                            error(message('stats:classreg:regr:LinearFormula:BadPredictorVarLength'));
                        end
                    else
                        error(message('stats:classreg:regr:LinearFormula:BadPredictorVarType'));
                    end
                    if haveVarNames
                        if isempty(responseName) % responseVar not given, or given as an index
                            if isscalar(respLoc)
                                responseName = varNames{respLoc};
                            elseif sum(predVars) == nvars-1
                                responseName = varNames{~predVars};
                            else
                                error(message('stats:classreg:regr:LinearFormula:AmbiguousResponse'));
                            end
                        else % responseVar given as a name
                            respLoc = find(strcmp(responseName,varNames));
                            if isempty(respLoc)
                                error(message('stats:classreg:regr:LinearFormula:BadResponseName'));
                            elseif predVars(respLoc)
                                error(message('stats:classreg:regr:LinearFormula:ResponseIsPredictor'));
                            end
                        end
                    else
                        if isempty(respLoc)
                            respLoc = find(~predVars);
                            if ~isscalar(respLoc)
                                error(message('stats:classreg:regr:LinearFormula:AmbiguousResponse'));
                            end
                        else
                            if predVars(respLoc)
                                error(message('stats:classreg:regr:LinearFormula:ResponseIsPredictor'));
                            end
                        end
                        if isempty(responseName)
                            responseName = 'y';
                        end
                        varNames = {};
                        varNames([1:(respLoc-1) (respLoc+1):nvars]) = internal.stats.numberedNames('x',1:nvars-1);
                        varNames{respLoc} = responseName;
                    end
                else % modelSpec == alias
                    alias = modelSpec;
                    if haveVarNames
                        if isempty(responseName) % responseVar not given, or given as an index
                            if isempty(respLoc)
                                respLoc = nvars;
                            end
                            responseName = varNames{respLoc};
                        else % responseVar given as a name
                            respLoc = find(strcmp(responseName,varNames));
                            if isempty(respLoc)
                                error(message('stats:classreg:regr:LinearFormula:BadResponseName'));
                            end
                        end
                    else
                        if isempty(responseName)
                            responseName = 'y';
                        end
                        if isempty(respLoc)
                            respLoc = nvars;
                        end
                        varNames = {};
                        varNames([1:(respLoc-1) (respLoc+1):nvars]) = internal.stats.numberedNames('x',1:nvars-1);
                        varNames{respLoc} = responseName;
                    end
                    predVars = true(1,nvars); predVars(respLoc) = false;
                end
                if isMissingArg(hasIntercept)
                    hasIntercept = true;
                elseif ~islogical(hasIntercept) || ~isscalar(hasIntercept)
                    error(message('stats:classreg:regr:LinearFormula:BadHasIntercept'));
                end
                if isMissingArg(link)
                    link = 'identity';
                end
                f.terms = classreg.regr.modelutils.model2terms(alias,predVars,hasIntercept);
                f.responseName = responseName;
                f.varNames = varNames;
                f.link = link;
                f.FunctionCalls = cell(1,0);
                f.str = composeFormulaString(f.link,f.responseName,terms2expr(f.terms,f.varNames));
                
            elseif internal.stats.isString(modelSpec)
                % The variables (both predictor and response) that occur in a
                % formula string define the default universe of variables, but
                % the universe may optionally be expanded with varNames.  Not
                % all variables in varNames need occur in the formula.
                [f.str,f.tree] = parseStr(f,modelSpec);
                if f.tree(1,1) ~= f.irules.Formula
                     error(message('stats:classreg:regr:LinearFormula:InvalidFormula', modelSpec));
                end

                if isMissingArg(varNames)
                    [f.str,f.tree] = substituteForAllElse(f,{});
                    varNames = getPredictorAndResponseNames(f);
                elseif internal.stats.isStrings(varNames,true)
                    varNames = mustBeUnique(asCellStr(varNames));
                    if ~isempty(setdiff(getPredictorAndResponseNames(f),varNames))
                        error(message('stats:classreg:regr:LinearFormula:BadFormulaVariables'));
                    end
                    [f.str,f.tree] = substituteForAllElse(f,varNames);
                else
                    error(message('stats:classreg:regr:LinearFormula:BadVarNameValue'));
                end
                f = processFormula(f,varNames);
                
                % These are inputs that ordinarily would not be specified,
                % because they're implicit in the model formula string.  Allow
                % them if they don't conflict with the string.
                if ~isMissingArg(responseVar)
                    if (internal.stats.isString(responseVar) && ~isequal(responseVar,f.responseName)) || ...
                       (internal.stats.isScalarInt(responseVar,1,length(varNames)) ...
                            && ~isequal(responseVar,find(strcmp(f.responseName,varNames))))
                        error(message('stats:classreg:regr:LinearFormula:ResponseVarConflict'));
                    end
                end
                if ~isMissingArg(hasIntercept) && ~isequal(hasIntercept,f.HasIntercept)
                    % Allow explicit request to override the formula
                    if f.HasIntercept
                        f = removeTerms(f,'1');
                    else
                        f = addTerms(f,'1');
                    end
                end
                if ~isMissingArg(link) && ~isequal(link,f.link)
                    if isequal(f.link,'identity')
                        % Link argument val used when no link in formula
                        f.link = link;
                        f.str = composeFormulaString(f.link,f.responseName,terms2expr(f.terms,f.varNames));
                    else
                        % Otherwise make sure the link and formula agree
                        error(message('stats:classreg:regr:LinearFormula:LinkConflict'));
                    end
                end
            else
                error(message('stats:classreg:regr:LinearFormula:BadModelSpec'));
            end
        end
        
        function f = addTerms(f,update)
            ntermsBefore = size(f.terms,1);
            if classreg.regr.FormulaProcessor.isTermsMatrix(update)
                f.terms = removeDupRows([f.terms; update]);
            elseif internal.stats.isString(update)
                [ustr,utree] = parseStr(f,update);
                if utree(1,1) ~= f.irules.LinearPredictor
                    error(message('stats:classreg:regr:LinearFormula:UpdateNotFormula'));
                end
                [addTerms,removeTerms] = createTerms(f,utree,ustr,1,f.varNames,false);
                f.terms = sortTerms(setdiff([f.terms; addTerms],removeTerms,'rows'));
            else
                error(message('stats:classreg:regr:LinearFormula:BadModelSpecUpdate'));
            end
            ntermsAfter = size(f.terms,1);
            if ntermsAfter<=ntermsBefore
                warning(message('stats:classreg:regr:LinearFormula:NoNewTerms'));
            end
            respLoc = find(strcmp(f.responseName,f.varNames));
            if ~isempty(respLoc) && any(f.terms(:,respLoc(1)))
                warning(message('stats:classreg:regr:LinearFormula:ResponseTerm'))
            end

            f.str = getStr(f);
        end
        function f = removeTerms(f,update)
            ntermsBefore = size(f.terms,1);
            if classreg.regr.FormulaProcessor.isTermsMatrix(update)
                [isfound,foundrow] = ismember(update,f.terms,'rows');
                f.terms(foundrow(isfound),:) = [];
            elseif internal.stats.isString(update)
                [ustr,utree] = parseStr(f,update);
                if utree(1,1) ~= f.irules.LinearPredictor
                    error(message('stats:classreg:regr:LinearFormula:UpdateNotFormula'));
                end
                [addTerms,removeTerms] = createTerms(f,utree,ustr,1,f.varNames,false);
                oldterms = [f.terms; removeTerms];
                f.terms = sortTerms(setdiff(oldterms,addTerms,'rows'));
            else
                error(message('stats:classreg:regr:LinearFormula:BadModelSpecUpdate'));
            end
            ntermsAfter = size(f.terms,1);
            if ntermsAfter>=ntermsBefore
                warning(message('stats:classreg:regr:LinearFormula:TermNotFound'));
            end
            f.str = getStr(f);
        end
        
        function fstr = char(f,maxWidth)
            fstr = prettyStr(f.str);
            if nargin == 2 && (length(fstr) > maxWidth)
                fstr = sprintf('%s',getString(message('stats:classreg:regr:LinearFormula:display_LinearFormula', ...
                    f.ResponseName,f.NTerms,f.NPredictors)));
% alternately ...
%                 start = regexp(fstr,' [+-] ','start');
%                 i = start(find(start < maxWidth,1,'last'));
%                 fstr = [fstr(1:i+2) '...'];
            end
        end
        
        function disp(f)
            fstr = prettyStr(f.str);
            strLHSlen = regexp(fstr,'\ *~\ *','start') - 1;
            
            % Fit the string to the command window width
            lf = sprintf('\n');
            pad = repmat(' ',1,strLHSlen+6);
            maxWidth = get(0,'CommandWindowSize'); maxWidth = maxWidth(1) - 1;
            start = regexp(fstr,' [+-] ','start');
            loc = 0;
            indent = 0;
            while length(fstr)-loc > maxWidth-indent
                i = start(find(start-loc < maxWidth-indent,1,'last'));
                fstr(i) = lf;
                loc = i + 1;
                indent = length(pad);
            end
            fstr = regexprep(fstr,'\n',[lf pad]);
            
            disp(fstr);
        end
        function f = removeCategoricalPowers(f,isCat,silent)
            quadCatTerms = find(max(f.Terms(:,isCat),[],2)>1);
            if ~isempty(quadCatTerms)
                if ~silent
                    warning(message('stats:classreg:regr:TermsRegression:NoCatPowers'));
                end
                f = removeTerms(f,f.Terms(quadCatTerms,:));
            end
        end
    end % public methods
    
    methods(Access='protected')
        function str = getStr(f)
            str = composeFormulaString(f.link,f.responseName,terms2expr(f.terms,f.varNames));
        end
        
        function [s,t] = substituteForAllElse(f,varNames)
            t = f.tree;
            s = f.str;

            % Substitute variable names for an "everything else" placeholder
            phLoc = find(t(1,:)==f.irules.AllElse);
            if isscalar(phLoc)
                if isMissingArg(varNames)
                    error(message('stats:classreg:regr:LinearFormula:BadAllElse',FormulaProcessor.allElseStr));
                end
                % "Everything else" is the set of variables that don't
                % explicitly appear in the formula
                explicitNames = getPredictorAndResponseNames(f);
                phStr = internal.stats.strCollapse(setdiff(varNames,explicitNames),' + ');
                if isempty(phStr)
                    phStr = '0';
                end
                modelString = [s(1:(t(2,phLoc)-1)) phStr s((t(3,phLoc)+1):end)];
                [s,t] = parseStr(f,modelString);
                if t(1,1) ~= f.irules.Formula
                     error(message('stats:classreg:regr:LinearFormula:InvalidFormula', modelString));
                end
            elseif ~isempty(phLoc)
                error(message('stats:classreg:regr:LinearFormula:MultipleAllElse', FormulaProcessor.allElseStr));
            end
        end
        
        function f = processFormula(f,varNames)
            t = f.tree;
            s = f.str;
            if t(1,2) ~= f.irules.Response
                error(message('stats:classreg:regr:LinearFormula:ResponseAndPredictor'));
            end

            % Identify the link function, if any
            if isfield(f.irules,'LinkFunction')
                j = find(t(1,:)==f.irules.LinkFunction,1,'first');
            else
                j = [];
            end
            if isempty(j)
                f.link = 'identity';
            else
                j = j+1; % get just the function name
                if t(1,j)==f.irules.PowerLink
                    k = j + find(t(1,(j+1):end)==f.irules.Expon,1,'first');
                    expon = s(t(2,k):t(3,k));
                    f.link = str2double(expon);
                    if isempty(f.link)
                        error(message('stats:classreg:regr:LinearFormula:UnrecognizedExponent', expon));
                    end
                else
                    f.link = s(t(2,j):t(3,j));
                end
            end

            % Identify the response
            j = find(t(1,:) == f.irules.ResponseVar);
            f.responseName = s(min(t(2,j)):max(t(3,j)));

            % Identify any function calls
%             j = find(tree(1,:) == f.irules.FunName);
%             names = cell(size(j));
%             for i = 1:length(j)
%                 names{i} = str(tree(2,j(i)):tree(3,j(i)));
%             end
%             f.FunctionCalls = uniqueLocal(names); % sorted

            % Identify terms in the linear predictor
            r = 2; % start node of the response subtree
            lp = r + t(5,r); % start node in LP subtree
            f.terms = sortTerms(createTerms(f,t,s,lp,varNames));
            f.varNames = varNames;
            f.str = getStr(f);
        end
        
        function names = getPredictorAndResponseNames(f)
            t = f.tree;
            j = find((t(1,:) == f.irules.ResponseVar) | ...
                     (t(1,:) == f.irules.PredictorVar));
%             j = find((t(1,:) == f.irules.ResponseVar) | ...
%                      (t(1,:) == f.irules.MatlabExpr) | ...
%                      (t(1,:) == f.irules.PredictorVar));
            names = cell(size(j));
            for i = 1:length(j)
                names{i} = f.str(t(2,j(i)):t(3,j(i)));
            end
            names = uniqueLocal(names); % sorted
        end
        
        function names = getPredictorNames(f)
            t = f.tree;
            j = find((t(1,:) == f.irules.PredictorVar));
%             j = find((t(1,:) == f.irules.PredictorVar) | ...
%                      (t(1,:) == f.irules.MatlabExpr));
            names = cell(size(j));
            for i = 1:length(j)
                names{i} = f.str(t(2,j(i)):t(3,j(i)));
            end
            names = uniqueLocal(names); % sorted
        end
    end % protected methods
    
    methods(Static)
        function tf = isTermsMatrix(x)
            tf = isnumeric(x) && ismatrix(x) && all(x(:)>=0) && all(x(:)==round(x(:)));
        end
        
        function tf = isModelAlias(model)
            if ~internal.stats.isString(model)
                tf = false;
            else
                switch lower(model)
                    case {'constant' 'linear' 'interactions' 'purequadratic' 'quadratic'}
                        tf = true;
                    otherwise
                        [c,startLoc,endLoc] = regexp(lower(model),'poly(\d*)','tokens');
                        tf = ( isscalar(c) && (startLoc == 1) && (endLoc == length(model)) );
                end
            end
        end
    end % static methods
end

function [str,tree] = parseStr(f,str)
    % Parse the formula string
    str = strtrim(str);
   
    global runTedsParser;
    if(isempty(runTedsParser))
        runTedsParser=false;
    end
    
    if(~runTedsParser) 
       treestruct = f.p.parse(str);
    else
       treestruct = classreg.regr.myPEGparser(str,...
           f.p.rulemap);
   end
    
    tree = treestruct.tree;
    if isempty(tree) || tree(3,1) < length(str)
        error(message('stats:classreg:regr:LinearFormula:BadString', str));
    end
end

function [addTerms,removeTerms] = createTerms(f,tree,str,start,vnames,addint)
    nvars = length(vnames);
    removeTerms = zeros(0,nvars);
    addTerms = identifyTerms(f,start+1);
    
    % Usually add an intercept.  It must be removed explicitly.
    if nargin<6 || addint
        addTerms = [zeros(1,nvars); addTerms];
    end
    
    if nargout < 2
        % Remove terms that were subtracted from those that were added, return
        % only net terms
        addTerms = setdiff(addTerms,removeTerms,'rows');
    else
        % Return added terms and subtracted terms separately
    end

    function [terms,node] = identifyTerms(f,node)
        switch tree(1,node)
        case f.irules.Sum
            last = node + tree(5,node) - 1;
            [terms,node] = identifyTerms(f,node+1);
            while node <= last
                terms1 = terms;
                type = tree(1,node);
                [terms2,node] = identifyTerms(f,node);
                switch type
                case f.irules.Addend
                    terms = [terms1; terms2];
                case f.irules.Subend
                    terms = terms1;
                end
            end
        case f.irules.Addend
            [terms,node] = identifyTerms(f,node+1);
        case f.irules.Subend
            terms = zeros(0,nvars);
            [rterms,node] = identifyTerms(f,node+1);
            removeTerms = [removeTerms; rterms];
        case f.irules.Product
            last = node + tree(5,node) - 1;
            [terms,node] = identifyTerms(f,node+1);
            while node <= last
                terms1 = terms;
                [terms2,node] = identifyTerms(f,node);
                n1 = size(terms1,1); n2 = size(terms2,1);
                i1 = repmat((1:n1)',1,n2); i2 = repmat(1:n2,n1,1);
                terms = [terms1; terms2; terms1(i1,:)+terms2(i2,:)];
            end
        case f.irules.Inside
            [terms1,node] = identifyTerms(f,node+1);
            [terms2,node] = identifyTerms(f,node);
            n1 = size(terms1,1); n2 = size(terms2,1);
            i1 = repmat((1:n1)',1,n2); i2 = repmat(1:n2,n1,1);
            terms = [terms1; terms1(i1,:)+terms2(i2,:)];
        case f.irules.Interaction
            last = node + tree(5,node) - 1;
            [terms,node] = identifyTerms(f,node+1);
            while node <= last
                terms1 = terms;
                [terms2,node] = identifyTerms(f,node);
                n1 = size(terms1,1); n2 = size(terms2,1);
                i1 = repmat((1:n1)',1,n2); i2 = repmat(1:n2,n1,1);
                terms = terms1(i1,:)+terms2(i2,:);
            end
        case f.irules.Power
            [terms2,node] = identifyTerms(f,node+1);
            expon = str2double(str(tree(2,node):tree(3,node)));
            node = node + 1;
            terms = terms2;
            for i = 2:expon
                terms1 = terms;
                n1 = size(terms1,1); n2 = size(terms2,1);
                i1 = repmat((1:n1)',1,n2); i2 = repmat(1:n2,n1,1);
                terms = [terms1; terms1(i1,:)+terms2(i2,:)];
            end
        case f.irules.PredictorVar
            name = str(tree(2,node):tree(3,node));
            ivar = find(strcmp(name,vnames));
            if isempty(ivar)
                error(message('stats:classreg:regr:LinearFormula:UnrecognizedVariable', name));
            end
            terms = zeros(1,nvars); terms(end,ivar) = 1;
            node = node + 1;
        case f.irules.Intercept
            terms = zeros(1,nvars);
            node = node + 1;
        case f.irules.Zero
            terms = zeros(0,nvars);
            node = node + 1;
        otherwise
            error(message('stats:classreg:regr:LinearFormula:UnrecognizedRule', tree( 1, node )));
        end
    end
end

function modelStr = composeFormulaString(link,responseName,expr)
    if isMissingArg(link) || isequal(link,'identity')
        modelStr = sprintf('%s ~ %s',responseName,expr);
    elseif isnumeric(link)
        modelStr = sprintf('power(%s,%d) ~ %s',responseName,expr);
    elseif internal.stats.isString(link)
        modelStr = sprintf('%s(%s) ~ %s',link,responseName,expr);
    else
        modelStr = sprintf('link(%s) ~ %s',responseName,expr);
    end
end

function expr = terms2expr(terms,varNames)
    if ~isempty(terms)
        termNames = classreg.regr.modelutils.terms2names(terms,varNames);
        termNames(strcmp(termNames,'(Intercept)')) = {'1'};
        keep = true(size(termNames));
        mainEffectTerms = find(sum(terms,2) == 1);
        mainEffectVars = find(sum(terms(mainEffectTerms,:),1));
        interactions = find((sum(terms,2) == 2) & (sum(terms>0,2) == 2))';
        for k = interactions
            ij = find(terms(k,:));
            [tf,loc] = ismember(ij,mainEffectVars);
            if all(tf)
                termNames{k} = internal.stats.strCollapse(varNames(ij),'*');
                keep(mainEffectTerms(loc)) = false;
%             elseif tf(1)
%                 termNames{k} = utils.strCollapse(varNames(ij),'/');
%                 keep(mainEffectTerms(loc(1))) = false;
%             elseif tf(2)
%                 termNames{k} = utils.strCollapse(varNames(ij([2 1])),'/');
%                 keep(mainEffectTerms(loc(2))) = false;
            end
        end
        expr = internal.stats.strCollapse(termNames(keep),' + ');
    elseif size(terms,1) > 0
        % An empty terms matrix with a term but no variables is an
        % intercept-only model.
        expr = '1';
    else % 0-by-0, or 0-by-nvars
        % An empty terms matrix with variables but no terms is a constant
        % no-intercept model.
        expr = '0';
    end
end

function terms = sortTerms(terms)
    nvars = size(terms,2);
    terms = unique(terms,'rows'); % remove duplicate terms
    terms = sortrows(terms,nvars:-1:1); % sort by var
    [~,ord] = sortrows([sum(terms,2) max(terms,[],2)]); terms = terms(ord,:); % then sort by term order
end

function x = removeDupRows(x)
    % remove duplicate rows, leaving x in the same order
    [~,i] = unique(x,'rows','first');
    if length(i) < size(x,1)
        x = x(sort(i),:);
    end
end


function str = prettyStr(str)
str(isspace(str)) = [];
tilde = find(str=='~',1);
if isempty(tilde)
    tilde = 0;
end
part1 = str(1:tilde-1);
part2 = str(tilde+1:end);
str = [part1,' ~ ',regexprep(part2,{'~' '+' '-'},{' ~ ' ' + ' ' - '})];
end

function b = uniqueLocal(a)
% unique would return a 0x1 given a 1x0, avoid that
b = unique(a); b = b(:)';
end

function a = mustBeUnique(a)
if length(a) ~= length(unique(a))
    error(message('stats:classreg:regr:LinearFormula:RepeatedVariables'));
end
a = a(:)';
end

function c = asCellStr(c)
if ~iscell(c), c = {c}; end
end

function tf = isMissingArg(x)
tf = isequal(x,[]);
end
