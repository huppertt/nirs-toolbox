classdef (Sealed = true) NonLinearFormula

%   Copyright 2011-2013 The MathWorks, Inc.

    properties(Constant,GetAccess='protected')
        rules = {
            '  Formula         = (Response "~")? Expression'
            '  Response        = Name'
            '+ Expression      = Test'
            '  Test            = Sum (("=="/"<="/">="/"<"/">"/"~=") Sum)?'
            '  Sum             = Product (("+"/"-") Product)*'
            '  Product         = Power (("*"/".*"/"/"/"./") Power)*'
            '  Power           = Term ("^" Term)?'
            '  Term            = "(" Expression ")" / SignedFactor / Number' 
            '  SignedFactor    = ("+"/"-")* Factor'
            '  Factor          = Function / PredictorOrCoef / MatlabExpr'   
            '  Function        = FunName "(" ArgList ")"'
            '  FunName         = Name'
            '  ArgList         = Expression ("," Expression)*'
            '  PredictorOrCoef = Name'
            '- Name            = [A-Za-z_] [A-Za-z0-9_]*'
            '  MatlabExpr      = "[" [^#x5B#x5D]+ "]"' % [^#x5B#x5D] => anything but "[" or "]"
            '  Number          = ("+"/"-")? [0-9]+ ("." [0-9]*)? / ("+"/"-")? ("." [0-9]+)?'
            };
        p = internal.stats.PEG(classreg.regr.NonLinearFormula.rules);
        irules = rulemap(classreg.regr.NonLinearFormula.p);
    end
    
    properties(GetAccess='protected',SetAccess='protected')
        fun = []; % empty if created from a string, function_handle if created from a function
        funType = ''; % '', 'simple', 'nested', 'anonymous' 'opaque-anonymous'
        tree = [];
        str = 'y ~ 0';
        rname = 'y'; % name of response variable
        cnames = cell(1,0); % names of model coefficients
        pnames = cell(1,0); % names of predictor variables
        vnames = {'y'}; % names of data variables
    end
    
    properties(Dependent,GetAccess='public',SetAccess='public')
        CoefficientNames
        VariableNames % names of all known variables, whether in model or not
        ResponseName
    end
    properties(Dependent,GetAccess='public',SetAccess='protected')
        Expression % the thing on the RHS of the '~'
        ModelFun
        InModel % which vars are in model as a predictor?
        
        Names % coef, variable, and constant names
        ExpressionNames % coef, predictor var, and constant names on the RHS of the '~'
        PredictorNames % names of predictors actually in model
        ConstantNames
        FunctionCalls
        NumCoefficients
        NVars
        NumPredictors
    end
    
    methods
        function rname = get.ResponseName(f)
            rname = f.rname;
        end
        function f = set.ResponseName(f,rname)
            f = parseStr(f,composeFormulaString(rname,f.Expression));
        end
        function n = get.NumCoefficients(f)
            n = length(f.cnames);
        end
        function cnames = get.CoefficientNames(f)
            cnames = f.cnames; % orginal order, not sorted
        end
        function f = set.CoefficientNames(f,cnames)
            f = setAllNames(f,cnames);
        end
        function n = get.NVars(f)
            n = length(f.vnames);
        end
        function n = get.NumPredictors(f)
            n = length(f.pnames);
        end
        function vnames = get.VariableNames(f)
            vnames = f.vnames; % original order, not sorted
        end
        function f = set.VariableNames(f,vnames)
            f = setAllNames(f,[],vnames);
        end
        function pnames = get.PredictorNames(f)
            pnames = f.pnames; % original order, not sorted
        end
        function names = get.ConstantNames(f)
            exprNames = f.ExpressionNames;
            names = setdiff(exprNames,[f.cnames f.vnames]); % leftovers are constants
        end
        function names = get.InModel(f)
            names = ismember(f.vnames,f.pnames);
        end
        
        function expr = get.Expression(f)
            t = f.tree;
            j = find(t(1,:) == f.irules.Expression);
            expr = prettyStr(f.str(t(2,j):t(3,j)));
        end
        function fun = get.ModelFun(f)
            if isempty(f.fun)
                expr = f.Expression;
                expr = regexprep(expr,{'(\.)?\*' '(\.)?/' '(\.)?\^'},{'.*' './' '.^'});
                betaNames = strcat({'b('},num2str((1:f.NumCoefficients)'),{')'});
                coefNames = f.CoefficientNames;
                expr = regexprepLocal(expr,strcat('\<',coefNames,'\>'),betaNames,strcmp('b',coefNames));
                Xnames = strcat({'X(:,'},num2str((1:f.NumPredictors)'),{')'});
                predNames = f.PredictorNames;
                expr = regexprepLocal(expr,strcat('\<',predNames,'\>'),Xnames,strcmp('X',predNames));
                fun = eval(['@(b,X) ' expr]);
            else
                fun = f.fun;
            end
        end
        
        function names = get.Names(f)
            if isOpaqueFunType(f.funType)
                % There is no string or function body to parse.
                names = cell(1,0);
            else
                t = f.tree;
                j = find(t(1,:) == f.irules.Response | t(1,:) == f.irules.PredictorOrCoef);
                names = cell(size(j));
                for i = 1:length(j)
                    names{i} = f.str(t(2,j(i)):t(3,j(i)));
                end
                names = uniqueLocal(names); % sorted
            end
        end
        function names = get.ExpressionNames(f)
            if isOpaqueFunType(f.funType)
                % There is no string or function body to parse.
                names = cell(1,0);
            else
                t = f.tree;
                j = find(t(1,:) == f.irules.PredictorOrCoef);
                names = cell(size(j));
                for i = 1:length(j)
                    names{i} = f.str(t(2,j(i)):t(3,j(i)));
                end
                names = uniqueLocal(names); % sorted
            end
        end
        function names = get.FunctionCalls(f)
            t = f.tree;
            j = find(t(1,:) == f.irules.FunName);
            names = cell(size(j));
            for i = 1:length(j)
                names{i} = f.str(t(2,j(i)):t(3,j(i)));
            end
            names = uniqueLocal(names); % sorted
        end
    end
        
    methods(Access='public')
        function f = NonLinearFormula(model,coefNames,predictorNames,responseName,varNames,ncoefs)
            if nargin < 1, return, end

            if nargin < 2, coefNames = []; end
            if nargin < 3, predictorNames = []; end
            if nargin < 4, responseName = []; end
            if nargin < 5, varNames = []; end
            if nargin < 6, ncoefs = 0; end
            
            givenFun = isa(model,'function_handle');
            haveAnonymous = false;
            if isa(model,'classreg.regr.NonLinearFormula')
                inferNames = false;
                f = model;
                % Formulas created from other formulas can modify the existing
                % response name.
                if ~isMissingArg(responseName)
                    f = parseStr(f,composeFormulaString(responseName,f.Expression));
                end
            else
                inferNames = true;
                fstr = '';
                if givenFun
                    funs = functions(model);
                    f.fun = model;
                    f.funType = funs.type;
                    
                    % If the model definition is a function, and response
                    % and/or predictor names were not given, and there are
                    % names for the data variables, go to those for response
                    % and predictor names.  If the model definition is a
                    % formula string or formula object, all names will be
                    % taken from the definition itself.
                    if ~isMissingArg(varNames)
                        if isMissingArg(responseName)
                            if isMissingArg(predictorNames)
                                responseName = varNames{end};
                                predictorNames = varNames(1:end-1);
                            else
                                responseName = setdiff(varNames,predictorNames);
                                if ~isscalar(responseName)
                                    error(message('stats:classreg:regr:NonLinearFormula:AmbiguousResponse'));
                                end
                                responseName = responseName{1};
                            end
                        elseif isMissingArg(predictorNames)
                            predictorNames = setdiff(varNames,responseName);
                        end
                    end
                    
                    % For non-anonymous functions, there's no function body to
                    % look at, and so no way to glean the coef and predictor
                    % names.  We could default to names like 'b1' and 'x1',
                    % but we don't even know how many there are.  So both sets
                    % of names must be provided for non-anonymous functions.
                    % Coef and predictor names can be gleaned from the arg
                    % list and body of an anonymous function, but also
                    % can be overridden.
                    if isOpaqueFunType(f.funType)
                        if isMissingArg(coefNames) || isMissingArg(predictorNames)
                            error(message('stats:classreg:regr:NonLinearFormula:PredictorCoefRequired', f.funType));
                        end
                    else
                        haveAnonymous = true;
                        if nargin(model)~=2
                            error(message('stats:classreg:regr:NonLinearFormula:BadModelArgs'));
                        end
                    end
                elseif internal.stats.isString(model)
                    fstr = model;
                else
                    error(message('stats:classreg:regr:NonLinearFormula:BadModelSpecification'));
                end
                try
                    if isempty(fstr)
                        % Get the formula expression from the function.  If it's
                        % an anonymous function, also get the coef and predictor
                        % names if they have not been provided explicitly.
                        [expr,coefNames,predictorNames] = function2expression(funs,coefNames,predictorNames);
                        if isMissingArg(responseName), responseName = 'y'; end
                        fstr = composeFormulaString(responseName,expr);
                    end
                    f = parseStr(f,fstr);
                catch ME
                    % If we ran into trouble deciphering an anonymous
                    % function, we can choose to treat it as opaque and
                    % display it in the style   y ~ F(b,x)   where F
                    % is intended to represent a generic function
                    if haveAnonymous
                        i1 = find(funs.function=='(',1);
                        i2 = find(funs.function==')',1);
                        funargs = funs.function(i1+1:i2-1);
                        f = parseStr(f,sprintf('%s ~ F(%s)',responseName,funargs));
                        f.funType = 'opaque-anonymous';
                        if isempty(coefNames)
                            coefNames = internal.stats.numberedNames('b',1:ncoefs)';
                            if ~isempty(varNames)
                                coefNames = genvarname(coefNames,varNames);
                            end
                        end
                    else
                        rethrow(ME)
                    end
                end
            end
            f = setAllNames(f,coefNames,varNames,predictorNames,responseName,inferNames);
        end
        
        % ----------------------------------------------------------------------
        function fstr = char(f,maxWidth)
            fstr = prettyStr(f.str);
            if nargin == 2 && (length(fstr) > maxWidth)
                fstr = sprintf('%s',getString(message('stats:classreg:regr:NonLinearFormula:display_NonlinearFormula', ...
                    f.ResponseName,f.NumCoefficients,f.NumPredictors)));
            end
        end
        
        % ----------------------------------------------------------------------
        function disp(f)
            fstr = prettyStr(f.str);
            strLHSlen = regexp(fstr,'\ *~\ *','start') - 1;
            
            % Fit the string to the command window width
            lf = sprintf('\n');
            pad = repmat(' ',1,strLHSlen+6);
            maxWidth = matlab.desktop.commandwindow.size; maxWidth = maxWidth(1) - 1;
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
        
        % ----------------------------------------------------------------------
        function print(f)
            f.p.pretty(struct('tree',f.tree,'string',f.str))
        end
    end % public methods
    
    methods(Access='protected')
        function f = parseStr(f,fstr)
            % Parse the formula string
            fstr(isspace(fstr)) = [];
            f.str = fstr;
            treestruct = f.p.parse(f.str);
            f.tree = treestruct.tree;
            t = f.tree;
            if t(3,1) < length(f.str)
                error(message('stats:classreg:regr:NonLinearFormula:InvalidFormula', fstr));
            elseif t(1,1) ~= f.irules.Formula
                error(message('stats:classreg:regr:NonLinearFormula:MissingFormulaContents'));
            end
        end
        
        % ----------------------------------------------------------------------
        function f = setAllNames(f,coefNames,varNames,predictorNames,responseName,infer)
            if nargin < 3, varNames = []; end
            if nargin < 4, predictorNames = []; end
            if nargin < 5, responseName = []; end
            if nargin < 6, infer = false; end
            
            givenFun = ~isempty(f.fun); % creating from a formula string
            
            t = f.tree;
            j = find(t(1,:) == f.irules.Response);
            f.rname = f.str(t(2,j):t(3,j));
            exprNames = f.ExpressionNames;
            
            % Formulas created from strings must respect the response name in
            % the string.  Formulas created from functions can have the
            % default response variable name ('y') overridden, but it will
            % already have been inserted in the formula string, so no need to
            % check it.
            if ~isMissingArg(responseName)
                if ~givenFun && ~isequal(responseName,f.rname)
                    error(message('stats:classreg:regr:NonLinearFormula:ResponseFormulaConflict'));
                end
            end
            
            % If we have a list of data variables, use that to determine which
            % names in the expression are predictors.
            if ~isMissingArg(varNames)
                f.vnames = checkNames(varNames,'Variable');
                if isMissingArg(predictorNames)
                    [~,i] = intersect(varNames,exprNames);
                    predictorNames = varNames(sort(i)); % preserve the order
                end
            end
            
            % Formulas created from an anonymous function can have the coef and
            % predictor variable names that appear in that functions overridden.
            % Formulas created from non-anonymous functions have had the coef
            % and predictor names specified explicitly.  Formulas created from
            % strings or other formulas can specify which is which, but must
            % respect the names that are in the expression.
            if ~givenFun
                if ~isMissingArg(predictorNames) && ~all(ismember(predictorNames,exprNames))
                    error(message('stats:classreg:regr:NonLinearFormula:MissingPredictors'));
                elseif ~isMissingArg(coefNames) && ~all(ismember(coefNames,exprNames))
                    error(message('stats:classreg:regr:NonLinearFormula:MissingCoefficients'));
                end
            end
            
            if isMissingArg(coefNames)
                if isMissingArg(predictorNames)
                    if infer
                        % Assume that coefs will look like b1, b2, etc., and
                        % everything else in the RHS expression is a
                        % predictor.
                        notCoefs = cellfun('isempty',regexp(exprNames,'\<b\d+\>'));
                        f.pnames = exprNames(notCoefs);
                        f.cnames = exprNames(~notCoefs);
                    else
                        % neither specified, leave existing names alone
                    end
                else
                    f.pnames = checkNames(predictorNames,'Predictor');
                    if infer
                        f.cnames = setdiff(exprNames,f.pnames); % leftovers are coefs
                    else
                        % Remove any overlap of predictors from coefs, but
                        % otherwise leave coefs alone.  Thus, any former
                        % predictors will become constants.
                        f.cnames = setdiff(f.cnames,f.pnames);
                    end
                end
            else
                f.cnames = checkNames(coefNames,'Coefficient');
                if isMissingArg(predictorNames)
                    if infer
                        f.pnames = setdiff(exprNames,f.cnames); % leftovers are predictors
                    else
                        % Remove any overlap of coefs from predictors, but
                        % otherwise leave predictors alone.  Thus, any former
                        % coefs will become constants.
                        f.pnames = setdiff(f.pnames,f.cnames);
                    end
                else % both coefs and predictors specified
                    f.pnames = checkNames(predictorNames,'Predictor');
                    if any(ismember(f.pnames,f.cnames))
                        error(message('stats:classreg:regr:NonLinearFormula:PredictorCoefficientConflict'));
                    end
                end
            end
            
            if isMissingArg(varNames)
                f.vnames = [f.pnames f.rname];
            else
                if ~all(ismember(f.pnames,f.vnames))
                    error(message('stats:classreg:regr:NonLinearFormula:PredictorNotVariable'));
                end
                if ~any(strcmp(f.rname,f.vnames))
                    error(message('stats:classreg:regr:NonLinearFormula:ResponseNotVariable'));
                end
            end
        end
    end % protected methods
    
    methods(Static,Access='public')
        function tf = isOpaqueFun(fun)
            if isa(fun,'function_handle')
                s = functions(fun);
                tf = isOpaqueFunType(s.type);
            else
                tf = false;
            end
        end
    end % static public methods
end
        
    
function tf = isOpaqueFunType(type)
tf = ~isempty(type) && ~strcmp(type,'anonymous');
end

function fstr = composeFormulaString(responseName,expressionStr)
    fstr = [responseName ' ~ ' expressionStr];
end

function [expr,coefNames,predictorNames] = function2expression(funs,coefNames,predictorNames)
    funType = funs.type;
    fstr = funs.function; fstr(isspace(fstr)) = [];

    if strcmp(funType,'anonymous')
        % Anonymous functions provide the names of coefficients and
        % predictor variables in their arg list.  By looking at the
        % subscripting in the function body, we can also tell the
        % numbers of each.
        %
        % Split the string into something like '@(b,X)' and the
        % expression itself.
        tokens = regexp(fstr,'@\(\s*([a-zA-Z_]\w*)\s*,\s*([a-zA-Z_]\w*)\s*\)\s*(.+)','tokens');
        coefBase = tokens{1}{1};
        predictorBase = tokens{1}{2};
        expr = tokens{1}{3};

        % Find occurrences of things like 'b(1)' and 'X(:,1)', where
        % the indices are integer literals.
        coefIndices = regexp(expr,['\<' coefBase '\>\((\d)\)'],'tokens');
        predictorIndices = regexp(expr,['\<' predictorBase '\>\(:\,(\d)\)'],'tokens');

        % Find occurrences of things like 'b' and 'X' not followed by
        % things like '(1)' and '(:,1)'.
        coefBadTokens = regexp(expr,['\<' coefBase '\>(?!\(\d\))'],'start');
        predictorBadTokens = regexp(expr,['\<' predictorBase '\>(?!\(:\,\d\))'],'tokens');
        if isempty(predictorIndices)
            predictorIndices = predictorBadTokens;
            predictorBadTokens = {};
        end

        % If unsuccessful, the anonymous function string isn't
        % something we can recognize.
        if isempty(coefIndices) || isempty(predictorIndices) ...
                || ~isempty(coefBadTokens) || ~isempty(predictorBadTokens)
            error(message('stats:classreg:regr:NonLinearFormula:UnrecognizedCoefficientPredictor', fstr));
        end
        
        % If successful, construct coef names like b1 and predictor names like
        % X1 if none are specified.  The names will be contiguously numbered
        % even if some elements/columns don't appear in the expression.
        coefIndices = cellstr(num2str((1:max(str2num(char([coefIndices{:}]'))))','%-d'));
        coefNamesFound = strcat(coefBase,coefIndices)';
        if isMissingArg(coefNames)
            coefNames = coefNamesFound;
        elseif length(coefNames) ~= length(coefNamesFound)
            error(message('stats:classreg:regr:NonLinearFormula:BadCoefficientNameLength'));
        end
        predictorIndices = cellstr(num2str((1:max(str2num(char([predictorIndices{:}]'))))','%-d'));
        if length(predictorIndices)>1
            predictorNamesFound = strcat(predictorBase,predictorIndices)';
        else
            predictorNamesFound = {predictorBase};
        end
        if isMissingArg(predictorNames)
            % No predictor name specified, so use the name in the expression
            predictorNames = predictorNamesFound;
        elseif length(predictorNames) == length(predictorNamesFound)
            % We have one name for each one the expression requires
        else
            % We're okay if we have more names than we need; omit the extra
            ok = ismember(predictorNamesFound,predictorNames);
            if all(ok)
                predictorNames = predictorNamesFound;
            elseif length(predictorIndices)>length(predictorNames)
                error(message('stats:classreg:regr:NonLinearFormula:UnrecognizedPredictors', length( predictorNames ), length( predictorIndices )));
            end
        end

        % Now replace things like 'b(1)' and 'X(:,1)' in the expression with
        % the specified or constructed coef and predictor names.
        coefPats = strcat(coefBase,'\(',coefIndices,'\)')';
        expr = regexprep(expr,coefPats,coefNames);
        predictorPats = strcat(predictorBase,'\(:,',predictorIndices,'\)')';
        expr = regexprep(expr,predictorPats,predictorNames(1:length(predictorPats)));

        % Replace elementwise operators with scalar operators.
        expr = regexprep(expr,{'\.\*' '\./' '\.\^'},{'*' '/' '^'});

    else
        % Non-anonymous functions are opaque to us, can't get any
        % info about names or numbers of coefficients or predictor
        % variables.
        expr = [fstr '(b,X)'];
    end
end

function str = prettyStr(str)
   % a stub
   str(isspace(str)) = [];
   str = regexprep(str,{'~' '+' '-'},{' ~ ' ' + ' ' - '});
   str = regexprep(str,{'(\.)?*' '(\.)?/' '(\.)?\^'},{'*' '/' '^'});
end

function a = checkNames(a,which)
if ~iscellstr(a)
    error(message('stats:classreg:regr:NonLinearFormula:InvalidNames', which));
end
% make sure the list of names has no duplicates
if length(a) ~= length(unique(a))
    error(message('stats:classreg:regr:NonLinearFormula:NamesNotUnique', which));
end
a = a(:)';
end

function b = uniqueLocal(a)
% unique would return a 0x1 given a 1x0, avoid that
b = unique(a); b = b(:)';
end

function string = regexprepLocal(string,pattern,replace,doFirst)
% regexprep does the replacements sequentially, and if one of the patterns is
% also one of the replacements, it might find and replace a string that it has
% just inserted from a previous pattern.  Avoid that by doing the specified
% pattern first.
i = [find(doFirst) find(~doFirst)];
string = regexprep(string,pattern(i),replace(i));
end

function tf = isMissingArg(x)
tf = isequal(x,[]);
end
