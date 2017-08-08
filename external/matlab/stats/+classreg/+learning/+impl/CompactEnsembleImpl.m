classdef CompactEnsembleImpl
%CompactEnsembleImpl Ensemble implementation.

%   Copyright 2010-2011 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=public,Hidden=true)
        Trained = {};        
        Combiner = [];
    end
    
    methods(Hidden)
        function this = CompactEnsembleImpl(trained,combiner)
            this.Trained = trained;
            this.Combiner = combiner;
        end
    end

    methods
        function [varargout] = predictorImportance(this,varargin)
            % Init
            trained = this.Trained;
            if isempty(trained)
                varargout{1} = [];
                if nargout>1
                    varargout{2} = [];
                end
                return;
            end
            T = length(trained);
            W = [];
            if ismember('LearnerWeights',properties(this.Combiner))
                W = this.Combiner.LearnerWeights;
            end
            if isempty(W)
                W = ones(T,1);
            end
            
            % Look at the 1st learner in the ensemble
            cmp = trained{1};
            if isa(cmp,'classreg.learning.classif.CompactClassifByBinaryRegr')
                cmp = cmp.CompactRegressionLearner;
            end
            
            % Is this an ensemble of learners of the same type?
            cls = class(cmp);
            if numel(trained)>1
                cfun = @(x) istype(x,cls);
                tf = cellfun(cfun,trained(2:end));
                if any(~tf)
                    error(message('stats:classreg:learning:impl:CompactEnsembleImpl:predictorImportance:NonUniformEnsemble'));
                end
            end
            
            % Is this an ensemble of decision trees?
            istree = false;
            if isa(cmp,'classreg.learning.classif.CompactClassificationTree') ...
                    || isa(cmp,'classreg.learning.regr.CompactRegressionTree')
                istree = true;
            end
            
            % Predictor importance is provided for ensembles of trees only
            if ~istree
                error(message('stats:classreg:learning:impl:CompactEnsembleImpl:predictorImportance:NonTreeLearner'));
            end
            
            % Is this an ensemble of trees with surrogate split info stored?
            dosurr = false;
            assoc = [];
            if nargout>1
                if istree && ~isempty(cmp.SurrCutFlip)
                    dosurr = true;
                else
                    if istree
                        warning(message('stats:classreg:learning:impl:CompactEnsembleImpl:predictorImportance:NoSurrInfo'));
                    else
                        warning(message('stats:classreg:learning:impl:CompactEnsembleImpl:predictorImportance:No2ndOutputArg'));
                    end
                end
            end

            % Get predictor importance and optionally predictor association
            % from the 1st learner
            imp = W(1)*predictorImportance(cmp);
            if dosurr
                assoc = W(1)*meanSurrVarAssoc(cmp);
            end
            
            % Loop over learners
            for t=2:T
                cmp = trained{t};
                if isa(cmp,'classreg.learning.classif.CompactClassifByBinaryRegr')
                    cmp = cmp.CompactRegressionLearner;
                end
                imp = imp + W(t)*predictorImportance(cmp);
                if dosurr
                    assoc = assoc + W(t)*meanSurrVarAssoc(cmp);
                end
            end
            
            % Return importance and association
            varargout{1} = imp/sum(W);
            if nargout>1
                varargout{2} = assoc/sum(W);
            end
        end
        
        function this = sortLearnersByWeight(this)
            [alpha,sorted] = sort(this.Combiner.LearnerWeights,'descend');
            this.Combiner = clone(this.Combiner,alpha);
            this.Trained = this.Trained(sorted);
        end
        
        function this = removeLearners(this,idx)
            if ~isnumeric(idx) || ~isvector(idx) ...
                    || any(idx<=0) || any(isnan(idx)) || any(isinf(idx)) ...
                    || max(idx)>numel(this.Trained)
                error(message('stats:classreg:learning:impl:CompactEnsembleImpl:removeLearners:BadIdx',...
                    numel( this.Trained )));
            end
            alpha = this.Combiner.LearnerWeights;
            alpha(idx) = [];
            this.Combiner = clone(this.Combiner,alpha);
            this.Trained(idx) = [];
        end
    end    
end

function tf = istype(cmpobj,expectedType)
tf = false;
if isa(cmpobj,'classreg.learning.classif.CompactClassifByBinaryRegr')
    cmpobj = cmpobj.CompactRegressionLearner;
end
cls = class(cmpobj);
if strcmp(cls,expectedType)
    tf = true;
end
end

