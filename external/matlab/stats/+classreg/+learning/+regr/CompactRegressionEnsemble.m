classdef CompactRegressionEnsemble < ...
        classreg.learning.regr.RegressionModel & classreg.learning.ensemble.CompactEnsemble
%CompactRegressionEnsemble Compact regression ensemble.
%   CompactRegressionEnsemble is a set of trained weak learner models.
%   It can predict ensemble response for new data by aggregating
%   predictions from its weak learners.
%
%   CompactRegressionEnsemble properties:
%       PredictorNames        - Names of predictors used for this ensemble.
%       CategoricalPredictors - Indices of categorical predictors.
%       ResponseName          - Name of the response variable.
%       ResponseTransform     - Transformation applied to predicted regression response.
%       NumTrained            - Number of trained learners in the ensemble.
%       Trained               - Trained learners.
%       TrainedWeights        - Learner weights.
%       CombineWeights        - Prescription for combining weighted learner predictions.
%       UsePredForLearner     - Use predictors for learners.
%
%   CompactRegressionEnsemble methods:
%       loss                  - Regression loss.
%       predict               - Predicted response of this model.
%       predictorImportance   - Importance of predictors for this model.
%       removeLearners        - Remove learners from this ensemble.
%
%   See also RegressionEnsemble.

%   Copyright 2010-2013 The MathWorks, Inc.


    methods(Access=protected)
        function this = CompactRegressionEnsemble(dataSummary,responseTransform,usepredforlearner)
            this = this@classreg.learning.regr.RegressionModel(dataSummary,responseTransform);
            this = this@classreg.learning.ensemble.CompactEnsemble(usepredforlearner);
        end
        
        function r = response(this,X,varargin)
            r = classreg.learning.ensemble.CompactEnsemble.aggregatePredict(...
                X,this.Impl.Combiner,this.Impl.Trained,[],[],NaN,...
                'usepredforlearner',this.UsePredForLearner,varargin{:});
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.regr.RegressionModel(this,s);
            s = propsForDisp@classreg.learning.ensemble.CompactEnsemble(this,s);
        end
    end
     
    % Concrete methods
    methods
        function yfit = predict(this,X,varargin)
        %PREDICT Predict response of the ensemble.
        %   YFIT=PREDICT(ENS,X) returns predicted response YFIT for regression
        %   ensemble ENS and matrix of predictors X. Data X must be a numeric
        %   matrix of size N-by-P, where P is the number of predictors used for
        %   training this model. YFIT is a vector of type double with size(X,1)
        %   elements.
        %
        %   YFIT=PREDICT(ENS,X,'PARAM1',val1,'PARAM2',val2,...) specifies
        %   optional parameter name/value pairs:
        %       'useobsforlearner' - Logical matrix of size N-by-NumTrained, where
        %                            N is the number of observations in X and
        %                            NumTrained is the number of weak learners.
        %                            This matrix specifies what learners in the
        %                            ensemble are used for what observations. By
        %                            default, all elements of this matrix are set
        %                            to true.
        %       'learners'         - Indices of weak learners in the ensemble
        %                            ranging from 1 to NumTrained. Only these
        %                            learners are used for making predictions. By
        %                            default, all learners are used.
        %
        %   See also RegressionEnsemble, CompactRegressionEnsemble.
        
            if isempty(X)
                yfit = predictEmptyX(this,X);
                return;
            end
            yfit = predict@classreg.learning.regr.RegressionModel(this,X,varargin{:});
        end
        
        function l = loss(this,X,Y,varargin)
        %LOSS Regression error.
        %   L=LOSS(ENS,X,Y) returns mean squared error for ensemble ENS computed
        %   using matrix of predictors X and observed response Y. Data X must be a
        %   numeric matrix of size N-by-P, where P is the number of predictors used
        %   for training this model. Y must be a vector of floating-point numbers
        %   with N elements.
        %
        %   L=LOSS(ENS,X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
        %   parameter name/value pairs:
        %       'lossfun'          - Function handle for loss, or string
        %                            representing a built-in loss function.
        %                            Available loss functions for regression:
        %                            'mse'. If you pass a function handle FUN, LOSS
        %                            calls it as shown below:
        %                               FUN(Y,Yfit,W)
        %                            where Y, Yfit and W are numeric vectors of
        %                            length N. Y is observed response, Yfit is
        %                            predicted response, and W is observation
        %                            weights. Default: 'mse'
        %       'weights'          - Vector of observation weights. By default the
        %                            weight of every observation is set to 1. The
        %                            length of this vector must be equal to the
        %                            number of rows in X.
        %       'useobsforlearner' - Logical matrix of size N-by-NumTrained, where
        %                            N is the number of observations in X and
        %                            NumTrained is the number of weak learners.
        %                            This matrix specifies what learners in the
        %                            ensemble are used for what observations. By
        %                            default, all elements of this matrix are set
        %                            to true.
        %       'learners'         - Indices of weak learners in the ensemble
        %                            ranging from 1 to NumTrained. Only these
        %                            learners are used for making predictions. By
        %                            default, all learners are used.
        %       'mode'             - 'ensemble' (default), 'individual' or
        %                            'cumulative'. If 'ensemble', this method
        %                            returns a scalar value for the full ensemble.
        %                            If 'individual', this method returns a vector
        %                            with one element per trained learner. If
        %                            'cumulative', this method returns a vector in
        %                            which element J is obtained by using learners
        %                            1:J from the input list of learners.
        %
        %   See also RegressionEnsemble, CompactRegressionEnsemble, predict.            
            
            N = size(X,1);
            args = {                  'lossfun' 'weights'};
            defs = {@classreg.learning.loss.mse ones(N,1)};
            [funloss,W,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            [X,Y,W] = prepareDataForLoss(this,X,Y,W);
            l = classreg.learning.ensemble.CompactEnsemble.aggregateLoss(...
                this.NTrained,X,Y,W,[],funloss,this.Impl.Combiner,...
                @classreg.learning.ensemble.CompactEnsemble.predictOneWithCache,...
                this.Impl.Trained,[],[],this.PrivResponseTransform,NaN,...
                'usepredforlearner',this.UsePredForLearner,extraArgs{:});
        end
        
        function [varargout] = predictorImportance(this,varargin)
        %PREDICTORIMPORTANCE Estimates of predictor importance.
        %   IMP=PREDICTORIMPORTANCE(ENS) computes estimates of predictor importance
        %   for ensemble ENS by summing these estimates over all weak learners in
        %   the ensemble. The returned vector IMP has one element for each input
        %   predictor in the data used to train this ensemble. A high value
        %   indicates that this predictor is important for this ensemble.
        %
        %   [IMP,MA]=PREDICTORIMPORTANCE(ENS) for ensembles of decision trees also
        %   returns a P-by-P matrix with predictive measures of association for P
        %   predictors. Element MA(I,J) is the predictive measure of association
        %   averaged over surrogate splits on predictor J for which predictor I is
        %   the optimal split predictor. PREDICTORIMPORTANCE averages this
        %   predictive measure of association over all trees in the ensemble.
        %
        %   See also RegressionEnsemble, CompactRegressionEnsemble,
        %   classreg.learning.classif.CompactRegressionTree/predictorImportance,
        %   classreg.learning.classif.CompactRegressionTree/meanSurrVarAssoc.
            
            [varargout{1:nargout}] = predictorImportance(this.Impl,varargin{:});
        end
    end

end
