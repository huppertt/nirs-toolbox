classdef (Sealed = true) MultivariateLinearFormula < classreg.regr.FormulaProcessor

%   Copyright 2011 The MathWorks, Inc.

    properties(Constant,GetAccess='protected')
        rules = {
            '  Formula          = (Response Spaces "~" Spaces)? LinearPredictor'
            '+ Response         = ResponseRange (Spaces "," Spaces ResponseRange)*'
            '  ResponseRange    = ResponseVar (Spaces "-" Spaces ResponseVar)?'
            '  ResponseVar      = Name'
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
%             '  FunName          = Name'
%             '  ArgList          = PredictorVar' % only simple elementwise calls
            '  PredictorVar     = Name'
            '- Name             = [A-Za-z_] [A-Za-z0-9_]*'
%             '  MatlabExpr       = "[" [^#x5B#x5D]+ "]"' % [^#x5B#x5D] => anything but "[" or "]"
            '  Expon            = [0-9]+ ( "." [0-9]+ )' %("+" / "-") [0-9]+ ( "." [0-9]+ )'
            '  Integer          = [0-9]+'
            '  Intercept        = "1"'
            '  Zero             = "0"'
           ['  AllElse          = "' classreg.regr.MultivariateLinearFormula.allElseStr '"']
            '- Spaces           = (" ")*'
            };
        p = internal.stats.PEG(classreg.regr.MultivariateLinearFormula.rules);
        irules = rulemap(classreg.regr.MultivariateLinearFormula.p);
    end
    properties(Access='protected')
        isMultivariate = true;
    end
    
    methods(Access='public')
        function f = MultivariateLinearFormula(varargin)
            f = f@classreg.regr.FormulaProcessor(varargin{:});
        end
    end
end
