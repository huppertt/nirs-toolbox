classdef (Sealed = true) LinearMixedFormula
%LinearMixedFormula A linear or a generalized linear mixed effects formula.
%   LinearMixedFormula represents a formula that defines a linear or a
%   generalized linear mixed effects regression model.
%
%   LinearMixedFormula methods:
%       LinearMixedFormula - constructor 
%       char - character array representation of LinearMixedFormula               
%       disp - display LinearMixedFormula 
%
%   LinearMixedFormula properties:
%       FELinearFormula - Information on the fixed effects part
%       GroupingVariableNames - Information on grouping variables
%       RELinearFormula - Information on random effects parts
%       ResponseName - Name of response
%       PredictorNames - Names of predictors in the formula
%       VariableNames - Names of variables used to create the formula
%       Link - Name of the link function 
%
%   See also classreg.regr.LinearFormula.
        
%   Copyright 2012-2013 The MathWorks, Inc.

    properties(Constant,GetAccess='protected')
         allElseStr = '...';
%rules - Rules defining the structure of a LinearMixedFormula string.
%   Compared to a LinearFormula, a LinearMixedFormula has two new rules:
%   LinearRandomPredictor and GroupingVar. A + sign on the left of both
%   these rules indicates that we want to record matches for both these
%   rules in the parse tree. The grouping variable can be a single
%   predictor variable or a multiway interaction between them. Grouping
%   variables cannot be continuous.
         rules = {
            '  Formula          = (Response Spaces "~" Spaces)? LinearPredictor Spaces ("+" LinearRandomPredictor)*'
            '+ Response         = LinkFunction / ResponseVar'
            '+ LinkFunction     = LinkName "(" Spaces ResponseVar Spaces ")" / PowerLink'
            '  PowerLink        = "power" "(" Spaces ResponseVar Spaces "," Spaces Expon Spaces ")"'
            '  ResponseVar      = Name'
            '+ LinearRandomPredictor = Spaces "(" Spaces LinearPredictor Spaces "|" Spaces GroupingVar Spaces ")" Spaces'
            '+ LinearPredictor  = Sum / Zero'
            '  Sum              = Augend (Spaces (Addend / Subend))*'
            '- Augend           = Conditional / Subend'
            '+ Addend           = "+" Spaces Conditional'
            '+ Subend           = "-" Spaces Conditional'
%             '  Conditional      = Interaction (Spaces "|" Spaces Interaction)?'
            '  Conditional      = Product'
            '  Product          = Inside (Spaces "*" Spaces Inside)*'
            '  Inside           = Interaction (Spaces "/" Spaces PredictorVar)?'
            '  Interaction      = Power (Spaces ":" Spaces Power)*'
            '  Power            = Predictor ("^" Integer)?'
%             '- Predictor        = "(" Spaces Sum Spaces ")" / AllElse / Function / PredictorVar / MatlabExpr / Intercept'
            '- Predictor        = "(" Spaces Sum Spaces ")" / AllElse / PredictorVar / Intercept'
%             '  Function         = FunName "(" Spaces ArgList Spaces ")"'
            '  LinkName         = Name'
%             '  FunName          = Name'
%             '  ArgList          = PredictorVar' % only simple elementwise calls
            '+  GroupingVar     = PredictorVar (Spaces ":" Spaces PredictorVar)*'
            '  PredictorVar     = Name'
            '- Name             = [A-Za-z_] [A-Za-z0-9_]*'
%             '  MatlabExpr       = "[" [^#x5B#x5D]+ "]"' % [^#x5B#x5D] => anything but "[" or "]"
            '  Expon            = [0-9]+ ( "." [0-9]+ )' %("+" / "-") [0-9]+ ( "." [0-9]+ )'
            '  Integer          = [0-9]+'
            '  Intercept        = "1"'
            '  Zero             = "0"'
           ['  AllElse          = "' classreg.regr.LinearMixedFormula.allElseStr '"']
            '- Spaces           = (" ")*'
            };
%p - A PEG parser for strings according to definitions given in rules.        
        p = internal.stats.PEG(classreg.regr.LinearMixedFormula.rules);
%irules - A structure mapping rule names to integers.        
        irules = rulemap(classreg.regr.LinearMixedFormula.p);
    end % protected constants
    
    properties(GetAccess='public',SetAccess='protected')
                
%FELinearFormula - A LinearFormula object representing the fixed effects
%   design matrix of the model.
        FELinearFormula
   
%GroupingVariableNames - A cell array containing grouping variable
%   information. GroupingVariableNames is of size G by 1 where G is the
%   number of random effects specifications in the model. Element i of
%   GroupingVariableNames is a cell array of length 1 by L such that the i
%   th grouping variable is an interaction between grouping variables
%   GroupingVariableNames{i}{1} through GroupingVariableNames{i}{L}.
        GroupingVariableNames
        
%RELinearFormula - A cell array of LinearFormula objects representing the
%   random effects design matrices of the model. If GroupingVariableNames
%   is of size G by 1 then RELinearFormula is also of size G by 1. Element
%   i of RELinearFormula contains a specification of the random effects
%   design matrix corresponding to grouping variable defined by elements of
%   GroupingVariableNames{i}.
        RELinearFormula
       
%ResponseName - A string containing the name of the response variable.        
        ResponseName

%PredictorNames - A cell array of strings containing the predictor names
%   for the response variable ResponseName.
        PredictorNames
        
%VariableNames - A cell array of strings containing the variable names
%   including the response variable. There may be some variables in
%   VariableNames that are not PredictorNames.
        VariableNames

%Link - A string containing the name of a link function. For linear mixed
%   effects models, Link will be equal to 'identity'.
        Link

    end % public properties
    
    properties(GetAccess='protected',SetAccess='protected')
%str - String representing the LinearMixedFormula object returned by the
%   parseStr method.
        str = 'y ~ 0';
%tree - Parse tree for str returned by the parseStr method.
        tree = [];
    end % protected properties

    methods(Access='public')
        
        function f = LinearMixedFormula(modelSpec,varNames)
%LinearMixedFormula Create a LinearMixedFormula object.
%   F = LinearMixedFormula(MODELSPEC) creates a LinearMixedFormula object F
%   from a string MODELSPEC representing a linear or a generalized linear
%   mixed effects model.
%
%   F = LinearMixedFormula(MODELSPEC,VARNAMES) also accepts a cell array of
%   strings VARNAMES. Variable names represented by strings in VARNAMES
%   must include those that appear in the MODELSPEC. All variable names
%   represented by VARNAMES are included in LinearMixedFormula even if some
%   are not used in MODELSPEC.
%
%   Examples:
%      % A linear mixed effects model.
%      modelSpec = 'y ~ 1 + x + (1 + x | Subject)';
%      F = classreg.regr.LinearMixedFormula(modelSpec);
%
%      % Create varNames with variable 'z' which is not used in modelSpec.
%      varNames = {'x','y','Subject','z'};
%      F = classreg.regr.LinearMixedFormula(modelSpec,varNames); % OK
%
%      % Create varNames without variable 'x'.
%      varNames = {'y','Subject','z'};
%      F = classreg.regr.LinearMixedFormula(modelSpec,varNames); % ERROR
%
%      % Another linear mixed effects model.
%      modelSpec = 'y ~ 1 + z1 + (1 + (z1+z3*g2):g5 | g6) + (1 |g7:g8)';
%      F = classreg.regr.LinearMixedFormula(modelSpec);
%
%      % Another linear mixed effects model.
%      modelSpec = 'mylink(y) ~ 1 + x1 + (1 + x2|g2) + (1 + x3|g3:g4) + (x4*x5 - x4 - 1 | g4:g3:g2)';
%      F = classreg.regr.LinearMixedFormula(modelSpec);

%   See also classreg.regr.LinearFormula.

            if nargin < 1, return, end
            if nargin < 2, varNames = []; end
                
            if internal.stats.isString(modelSpec)
                % The variables (both predictor and response) that occur in a
                % formula string define the default universe of variables, but
                % the universe may optionally be expanded with varNames.  Not
                % all variables in varNames need occur in the formula.
                [f.str,f.tree] = parseStr(modelSpec);
                if f.tree(1,1) ~= classreg.regr.LinearMixedFormula.irules.Formula
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
            else
                error(message('stats:classreg:regr:LinearFormula:BadModelSpec'));
            end
        end

        function fstr = char(f,maxWidth)
%char Returns a character array representing the LinearMixedFormula.
%   FSTR = char(F) returns a character array that represents the 
%   LinearMixedFormula F.
%
%   FSTR = char(F,MAXWIDTH) accepts an integer specifying the maximum width
%   of the formula string. If the formula string is wider than MAXWIDTH, a
%   short textual description of the LinearMixedFormula is displayed.
            
            fstr = getDisplayString(f);       
                        
            if nargin == 2 && (length(fstr) > maxWidth)
                fstr = sprintf('%s',['Linear Mixed Formula with ',num2str(length(f.PredictorNames)),' predictors.']);
            end
            
        end
        
        function disp(f)
%disp Displays the LinearMixedFormula.
%   disp(F) displays the LinearMixedFormula F on the screen.

            fstr = getDisplayString(f);    
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
                
    end % public methods
    
    methods(Access='protected')
        
        function [s,t] = substituteForAllElse(f,varNames)
            t = f.tree;
            s = f.str;

            % Substitute variable names for an "everything else" placeholder
            phLoc = find(t(1,:)==classreg.regr.LinearMixedFormula.irules.AllElse);
            if isscalar(phLoc)
                if isMissingArg(varNames)
                    error(message('stats:classreg:regr:LinearFormula:BadAllElse', classreg.regr.LinearMixedFormula.allElseStr));
                end
                % "Everything else" is the set of variables that don't
                % explicitly appear in the formula
                explicitNames = getPredictorAndResponseNames(f);
                phStr = internal.stats.strCollapse(setdiff(varNames,explicitNames),' + ');
                if isempty(phStr)
                    phStr = '0';
                end
                modelString = [s(1:(t(2,phLoc)-1)) phStr s((t(3,phLoc)+1):end)];
                [s,t] = parseStr(modelString);
                if t(1,1) ~= classreg.regr.LinearMixedFormula.irules.Formula
                     error(message('stats:classreg:regr:LinearFormula:InvalidFormula', modelString));
                end
            elseif ~isempty(phLoc)
                error(message('stats:classreg:regr:LinearFormula:MultipleAllElse', classreg.regr.LinearMixedFormula.allElseStr));
            end
        end
        
        function f = processFormula(f,varNames)
            t = f.tree;
            s = f.str;
            if t(1,2) ~= classreg.regr.LinearMixedFormula.irules.Response
                error(message('stats:classreg:regr:LinearFormula:ResponseAndPredictor'));
            end
            
            %*************************************************************
            % Meaning of s and t:
            %
            % treestruct = classreg.regr.LinearMixedFormula.p.parse(f.str); 
            %
            %  s = treestruct.string = Formula string entered by the user
            %  t = treestruct.tree   = parse tree of recorded rules
            %
            %  t is a 5-by-N array of integers. Each column represents a 
            %  node in the tree with rows containing the following info:
            %      rule number for node
            %      start index of match for node
            %      end index of match for node
            %      offset to parent (subtract) of node
            %      number of nodes in the subtree rooted at the node
            %  Any child nodes start directly to the right of the parent
            %  node and appear in order left to right. Thus the entire
            %  subtree rooted at node N is between columns N and
            %  N+TREE(5,N)-1, inclusive.
            %**************************************************************
            
            % (1) Create a string representing the fixed effects spec and 
            %     then convert it into a LinearFormula.
            
                % Node number for rule Response
                idx = find(t(1,:) == classreg.regr.LinearMixedFormula.irules.Response,1,'first');
                % Part of s that matches Response
                ResponseStr = s(t(2,idx):t(3,idx));      
                % Append 'space ~ space' to ResponseStr
                FELinearFormulaStr = [ResponseStr,' ~ '];
                % Node number for rule LinearPredictor (for fixed effects)
                idx = find(t(1,:) == classreg.regr.LinearMixedFormula.irules.LinearPredictor,1,'first');
                % Append string representing LinearPredictor
                FELinearFormulaStr = [FELinearFormulaStr,s(t(2,idx):t(3,idx))]; 
                % Make f.FELinearFormula
                f.FELinearFormula = classreg.regr.LinearFormula(FELinearFormulaStr,varNames);
    
            % (2) Create a cell array of grouping variables and the 
            %     corresponding random effects formula strings for each 
            %     random effects spec.
            
                % Node numbers for LinearRandomPredictor
                lrp = find(t(1,:) == classreg.regr.LinearMixedFormula.irules.LinearRandomPredictor);
                % Node numbers for LinearPredictor
                lp = find(t(1,:) == classreg.regr.LinearMixedFormula.irules.LinearPredictor); 
                % Node numbers for GroupingVar
                gv = find(t(1,:) == classreg.regr.LinearMixedFormula.irules.GroupingVar);
                % Node numbers for PredictorVar
                pv = find(t(1,:) == classreg.regr.LinearMixedFormula.irules.PredictorVar);
                % Initialize RELinearFormulaStr and f.GroupingVariableNames
                RELinearFormulaStr = cell(length(lrp),1);
                f.GroupingVariableNames = cell(length(lrp),1);
                for i = 1:length(lrp)
                    % List of nodes that are children of this LinearRandomPredictor
                    childidx = lrp(i) : lrp(i) + t(5,lrp(i)) - 1;
                    % Node number for LinearPredictor inside LinearRandomPredictor
                    lpidx = intersect(childidx, lp);
                    % Append LinearPredictor for ith random effects spec
                    RELinearFormulaStr{i} = [ResponseStr,' ~ ',s(t(2,lpidx):t(3,lpidx))];
                    % Make f.RELinearFormula{i}
                    f.RELinearFormula{i} = classreg.regr.LinearFormula(RELinearFormulaStr{i},varNames);
                    % Node number for GroupingVar inside LinearRandomPredictor
                    gvidx = intersect(childidx, gv);                    
                        % List of nodes that are children of GroupingVar node.
                        gvchildidx = gvidx : gvidx + t(5,gvidx) - 1;
                            % Extract PredictorVars that make up the GroupingVar
                            pidx = intersect(gvchildidx, pv);
                            f.GroupingVariableNames{i} = cell(1,length(pidx));
                            for j = 1:length(pidx)
                                f.GroupingVariableNames{i}{j} = s(t(2,pidx(j)):t(3,pidx(j)));
                            end
                end
                
            % (3) Identify the link function, if any
            j = find(t(1,:)==classreg.regr.LinearMixedFormula.irules.LinkFunction,1,'first');
            if isempty(j)
                f.Link = 'identity';
            else
                j = j+1; % get just the function name
                if t(1,j)==classreg.regr.LinearMixedFormula.irules.PowerLink
                    k = j + find(t(1,(j+1):end)==classreg.regr.LinearMixedFormula.irules.Expon,1,'first');
                    expon = s(t(2,k):t(3,k));
                    f.Link = str2double(expon);
                    if isempty(f.link)
                        error(message('stats:classreg:regr:LinearFormula:UnrecognizedExponent', expon));
                    end
                else
                    f.Link = s(t(2,j):t(3,j));
                end
            end

            % (4) Identify the response
            j = find(t(1,:) == classreg.regr.LinearMixedFormula.irules.ResponseVar);
            f.ResponseName = s(t(2,j):t(3,j));

            % (5) Identify predictor names
            f.PredictorNames = getPredictorNames(f);
            
            % (6) Identify variable names
            f.VariableNames = varNames;
                        
        end
        
        function names = getPredictorAndResponseNames(f)
            t = f.tree;
            j = find((t(1,:) == classreg.regr.LinearMixedFormula.irules.ResponseVar) | ...
                     (t(1,:) == classreg.regr.LinearMixedFormula.irules.PredictorVar));
%             j = find((t(1,:) == classreg.regr.LinearMixedFormula.irules.ResponseVar) | ...
%                      (t(1,:) == classreg.regr.LinearMixedFormula.irules.MatlabExpr) | ...
%                      (t(1,:) == classreg.regr.LinearMixedFormula.irules.PredictorVar));
            names = cell(size(j));
            for i = 1:length(j)
                names{i} = f.str(t(2,j(i)):t(3,j(i)));
            end
            names = uniqueLocal(names); % sorted
        end
        
        function names = getPredictorNames(f)
            t = f.tree;
            j = find((t(1,:) == classreg.regr.LinearMixedFormula.irules.PredictorVar));
%             j = find((t(1,:) == classreg.regr.LinearMixedFormula.irules.PredictorVar) | ...
%                      (t(1,:) == classreg.regr.LinearMixedFormula.irules.MatlabExpr));
            names = cell(size(j));
            for i = 1:length(j)
                names{i} = f.str(t(2,j(i)):t(3,j(i)));
            end
            names = uniqueLocal(names); % sorted
        end
        
        function str = getDisplayString(f)            
            str =  char(f.FELinearFormula);
            gnames = prettyGroupingVariableNames(f.GroupingVariableNames);
             for i = 1:length(f.RELinearFormula)
                appendStr = ['(',f.RELinearFormula{i}.LinearPredictor,' | ',gnames{i},')'];
                str = [str,' + ',appendStr]; %#ok<AGROW>
             end            
        end
                
    end % protected methods
    
end

function prettynames = prettyGroupingVariableNames(names)
% names is a G by 1 cell array.
% names{i} is a 1 by L cell array. 
% If names{i} contains {'g1','g2'} then prettynames{i} will contain 
% 'g1:g2'. prettynames will be a G by 1 cell array of strings.
    G = length(names);
    prettynames = cell(G,1);
    for i = 1:G
        prettynames{i} = names{i}{1};
        for j = 2:length(names{i})
           prettynames{i} = [prettynames{i},':',names{i}{j}]; 
        end        
    end
end

function [str,tree] = parseStr(str)
    % Parse the formula string
    str = strtrim(str);
   
   global runTedsParser;
    if(isempty(runTedsParser))
        runTedsParser=false;
    end
    
    if(~runTedsParser) 
       treestruct = classreg.regr.LinearMixedFormula.p.parse(str);
    else
        treestruct = classreg.regr.myPEGparser(str,...
           classreg.regr.LinearMixedFormula.p.rulemap);
   end

    
    tree = treestruct.tree;
     if isempty(tree) || tree(3,1) < length(str)
         error(message('stats:classreg:regr:LinearFormula:BadString', str));
     end
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
