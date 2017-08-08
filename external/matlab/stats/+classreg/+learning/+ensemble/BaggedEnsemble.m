classdef BaggedEnsemble
%BaggedEnsemble Ensemble grown by resampling.
%   BaggedEnsemble is the super class for ensemble models grown by
%   resampling the training data.
    
%   Copyright 2010-2013 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=protected,Abstract=true)
        ModelParams;
    end
    
    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        %FRESAMPLE Fraction of training data for resampling.
        %   The FResample property is a numeric scalar between 0 and 1. It is set
        %   to the fraction of the training data resampled at random for every weak
        %   learner in this ensemble.
        %
        %   See also classreg.learning.ensemble.BaggedEnsemble.
        FResample;
        
        %REPLACE Flag indicating if training data were sampled with replacement.
        %   The Replace property is a logical flag. It is set to true if the
        %   training data for weak learners in this ensemble were sampled with
        %   replacement and set to false otherwise.
        %
        %   See also classreg.learning.ensemble.BaggedEnsemble.
        Replace;
        
        %USEOBSFORLEARNER Use observations for learners.
        %   The UseObsForLearner property is a logical matrix of size
        %   N-by-NumTrained, where N is the number of observations in the training
        %   data and NumTrained is the number of trained weak learners. An element
        %   (I,J) of this matrix is set to true if observation I was used for
        %   training learner J and set to false otherwise.
        %
        %   See also classreg.learning.ensemble.BaggedEnsemble.
        UseObsForLearner;
    end
    
    methods(Access=protected)
        function this = BaggedEnsemble()
        end
        
        function s = propsForDisp(this,s)
            if nargin<2 || isempty(s)
                s = struct;
            else
                if ~isstruct(s)
                    error(message('stats:classreg:learning:ensemble:BaggedEnsemble:propsForDisp:BadS'));
                end
            end
            s.FResample = this.FResample;
            s.Replace = this.Replace;
            s.UseObsForLearner = this.UseObsForLearner;
        end
    end
    
    methods
        function fresample = get.FResample(this)
            fresample = this.ModelParams.Generator.FResample;
        end
        
        function replace = get.Replace(this)
            replace = this.ModelParams.Generator.Replace;
        end
        
        function usenfort = get.UseObsForLearner(this)
            usenfort = this.ModelParams.Generator.UseObsForIter;
        end
    end
    
end
