classdef ByBinaryRegrParams < classreg.learning.modelparams.ModelParams
%ByBinaryRegrParams Parameters for binary classification by regression.
%
%   ByBinaryRegrParams properties:
%       RegressionTemplate  - Learner template for regression.
    
%   Copyright 2010 The MathWorks, Inc.


    properties
        RegressionTemplate = [];
    end
    
    methods(Access=protected)
        function this = ByBinaryRegrParams(regtmp)
            this = this@classreg.learning.modelparams.ModelParams('ByBinaryRegr','classification');
            this.RegressionTemplate = regtmp;
        end
    end

    methods(Static,Hidden)
        function [holder,extraArgs] = make(type,varargin)
            % Decode input args
            args = {'learner'};
            defs = {       []};
            [regtmp,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Check input args
            if ~isempty(regtmp) && ~isa(regtmp,'classreg.learning.FitTemplate')
                error(message('stats:classreg:learning:modelparams:ByBinaryRegrParams:make:LearnerNotFitTemplate'));
            end
            if ~isempty(regtmp) && ~strcmp(regtmp.Type,'regression')
                error(message('stats:classreg:learning:modelparams:ByBinaryRegrParams:make:LearnerNotForRegression'));
            end

            % Make holder
            holder = classreg.learning.modelparams.ByBinaryRegrParams(regtmp);
        end
    end

    methods(Hidden)
        function this = fillDefaultParams(this,X,Y,W,dataSummary,classSummary)
            if isempty(this.RegressionTemplate)
                this.RegressionTemplate = classreg.learning.FitTemplate.make(...
                    'Tree','type','regression','minleaf',5);
            end
        end
    end

end
