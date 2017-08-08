classdef FullRegressionModel < ...
        classreg.learning.FullClassificationRegressionModel & classreg.learning.regr.RegressionModel
%FullRegressionModel Full regression model.
%   FullRegressionModel is the super class for full regression
%   models represented by objects storing the training data. This class is
%   derived from RegressionModel.
%
%   See also classreg.learning.regr.RegressionModel.

%   Copyright 2010-2014 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        %Y Observed response used to train this model.
        %   The Y property is a vector of type double.
        %
        %   See also classreg.learning.regr.FullRegressionModel.
        Y;
    end
    
    methods
        function y = get.Y(this)
            y = this.PrivY;
        end
    end
        
    methods(Access=protected)
        function this = FullRegressionModel(X,Y,W,modelParams,dataSummary,responseTransform)
            this = this@classreg.learning.FullClassificationRegressionModel(...
                dataSummary,X,Y,W,modelParams);
            this = this@classreg.learning.regr.RegressionModel(dataSummary,responseTransform);
            this.ModelParams = fillIfNeeded(modelParams,X,Y,W,dataSummary,[]);
        end
        
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.regr.RegressionModel(this,s);
            s = propsForDisp@classreg.learning.FullClassificationRegressionModel(this,s);
        end
    end
    
    methods
        function partModel = crossval(this,varargin)
        %CROSSVAL Cross-validate this model.
        %   CVMODEL=CROSSVAL(MODEL) builds a partitioned model CVMODEL from model
        %   MODEL represented by a full object for regression. You can then
        %   assess the predictive performance of this model on cross-validated data
        %   using methods and properties of CVMODEL. By default, CVMODEL is built
        %   using 10-fold cross-validation on the training data.
        %
        %   CVMODEL=CROSSVAL(MODEL,'PARAM1',val1,'PARAM2',val2,...) specifies
        %   optional parameter name/value pairs:
        %      'KFold'      - Number of folds for cross-validation, a numeric
        %                     positive scalar; 10 by default.
        %      'Holdout'    - Holdout validation uses the specified
        %                     fraction of the data for test, and uses the rest of
        %                     the data for training. Specify a numeric scalar
        %                     between 0 and 1.
        %      'Leaveout'   - If 'on', use leave-one-out cross-validation.
        %      'CVPartition' - An object of class CVPARTITION; empty by default. If
        %                      a CVPARTITION object is supplied, it is used for
        %                      splitting the data into subsets.
        %
        %   See also classreg.learning.regr.FullRegressionModel,
        %   cvpartition,
        %   classreg.learning.partition.RegressionPartitionedModel.

            idxBaseArg = find(ismember(varargin(1:2:end),...
                classreg.learning.FitTemplate.AllowedBaseFitObjectArgs));
            if ~isempty(idxBaseArg)
                error(message('stats:classreg:learning:regr:FullRegressionModel:crossval:NoBaseArgs', varargin{ 2*idxBaseArg - 1 }));
            end
            temp = classreg.learning.FitTemplate.make(this.ModelParams.Method,...
                'type','regression','responsetransform',this.PrivResponseTransform,...
                'modelparams',this.ModelParams,'CrossVal','on',varargin{:});
            partModel = fit(temp,this.X,this.Y,'Weights',this.W,...
                'predictornames',this.PredictorNames,'categoricalpredictors',this.CategoricalPredictors,...
                'responsename',this.ResponseName);
        end
        
        function [varargout] = resubPredict(this,varargin)
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            [varargout{1:nargout}] = predict(this,this.X,varargin{:});
        end
        
        function [varargout] = resubLoss(this,varargin)
            classreg.learning.FullClassificationRegressionModel.catchWeights(varargin{:});
            [varargout{1:nargout}] = ...
                loss(this,this.X,this.Y,'Weights',this.W,varargin{:});
        end        
    end

    methods(Static,Hidden)
        function [X,Y,W,dataSummary,responseTransform] = prepareData(X,Y,varargin)
            % Process input args
            args = {'responsetransform'};
            defs = {                 []};
            [transformer,~,crArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Pre-process
            if ~isfloat(X)
                error(message('stats:classreg:learning:regr:FullRegressionModel:prepareData:BadXType'));
            end
            [X,Y,W,dataSummary] = ...
                classreg.learning.FullClassificationRegressionModel.prepareDataCR(X,Y,crArgs{:});

            % Check Y type
            if ~isfloat(Y) || ~isvector(Y)
                error(message('stats:classreg:learning:regr:FullRegressionModel:prepareData:BadYType'));
            end
            internal.stats.checkSupportedNumeric('Y',Y,true);
            Y = Y(:);
            
            % Get rid of NaN's in response
            t = isnan(Y);
            if any(t)
                Y(t) = [];
                X(t,:) = [];
                W(t) = [];
            end
            if isempty(X)
                error(message('stats:classreg:learning:regr:FullRegressionModel:prepareData:NoGoodYData'));
            end
            
            % Renormalize weights
            W = W/sum(W);

            % Make output response transformation
            if isempty(transformer)
                responseTransform = @classreg.learning.transform.identity;
            elseif ischar(transformer)
                if strcmpi(transformer,'none')
                    responseTransform = @classreg.learning.transform.identity;
                else
                    responseTransform = str2func(['classreg.learning.transform.' transformer(:)']);
                end
            else
                if ~isa(transformer,'function_handle')
                    error(message('stats:classreg:learning:regr:FullRegressionModel:prepareData:BadResponseTransformation'));
                end
                responseTransform = transformer;
            end
         end
    end

end
