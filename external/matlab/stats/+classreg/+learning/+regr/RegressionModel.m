classdef RegressionModel < classreg.learning.Predictor

%   Copyright 2010-2014 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=protected,Hidden=true)
        PrivResponseTransform = [];
    end
    
    properties(GetAccess=public,SetAccess=public,Dependent=true)
        %RESPONSETRANSFORM Transformation applied to predicted response.
        %   The ResponseTransform property is a string describing how raw response
        %   values predicted by the model are transformed. You can assign a
        %   function handle to this property.
        %
        %   See also classreg.learning.regr.RegressionModel.
        ResponseTransform;
    end
    
    methods
        function ts = get.ResponseTransform(this)
            ts = func2str(this.PrivResponseTransform);
            if strcmp(ts,'classreg.learning.transform.identity')
                ts = 'none';
            end
            idx = strfind(ts,'classreg.learning.transform.');
            if ~isempty(idx)
               ts = ts(1+length('classreg.learning.transform.'):end); 
            end
        end
        
        function this = set.ResponseTransform(this,rt)
            this.PrivResponseTransform = ...
                classreg.learning.internal.convertScoreTransform(rt,'handle',1);
        end
    end
    
    methods(Hidden)
        function this = RegressionModel(dataSummary,responseTransform)
            this = this@classreg.learning.Predictor(dataSummary);
            this.PrivResponseTransform = responseTransform;
        end
    end
    
    methods(Access=protected)
        function s = propsForDisp(this,s)
            s = propsForDisp@classreg.learning.Predictor(this,s);
            s.ResponseTransform = this.ResponseTransform;
        end
        
        function yfit = predictEmptyX(this,X)
            D = numel(this.PredictorNames);
            if size(X,2)~=D
                error(message('stats:classreg:learning:regr:RegressionModel:predictEmptyX:XSizeMismatch', D));
            end
            yfit = NaN(0,1);
        end
        
        function [X,Y,W] = prepareDataForLoss(this,X,Y,W)
            % Check X
            if ~isnumeric(X) || ~ismatrix(X)
                error(message('stats:classreg:learning:regr:RegressionModel:prepareDataForLoss:BadXType'));
            end
            
            % Check Y type
            if ~isempty(Y) && (~isfloat(Y) || ~isvector(Y))
                error(message('stats:classreg:learning:regr:RegressionModel:prepareDataForLoss:BadYType'));
            end
            internal.stats.checkSupportedNumeric('Y',Y);
            Y = Y(:);

            % Check size
            N = size(X,1);
            if numel(Y)~=N
                error(message('stats:classreg:learning:regr:RegressionModel:prepareData:SizeXYMismatch'));
            end
            
            % Check weights
            if ~isfloat(W) || ~isvector(W) || length(W)~=N || any(W<0)
                error(message('stats:classreg:learning:regr:RegressionModel:prepareData:BadWeights', N));
            end
            internal.stats.checkSupportedNumeric('Weights',W,true);
            W = W(:);
            
            % Remove observations for NaN responses
            t = isnan(Y);
            if any(t) && N>0
                X(t,:) = [];
                Y(t,:) = [];
                W(t,:) = [];
            end
            
            % Normalize weights
            if sum(W)>0
                W = W/sum(W);
            end
        end
    end
       
    methods(Access=protected,Abstract=true)
        r = response(this,X,varargin)
    end
    
    methods
        function Yfit = predict(this,X,varargin)
        %PREDICT Predict response of the model.
        %   YFIT=PREDICT(MODEL,X) returns predicted response YFIT for regression
        %   model MODEL and matrix of predictors X. Data X must be a numeric
        %   matrix of size N-by-P, where P is the number of predictors used for
        %   training this model. YFIT is a vector of type double with size(X,1)
        %   elements.
        %
        %   See also classreg.learning.regr.RegressionModel,
        %   classreg.learning.regr.RegressionModel/predict.
            
            if isempty(X)
                Yfit = predictEmptyX(this,X);
                return;
            end
            Yfit = this.PrivResponseTransform(response(this,X,varargin{:}));
        end
        
        function l = loss(this,X,Y,varargin)
        %LOSS Regression error.
        %   L=LOSS(MODEL,X,Y) returns mean squared error for MODEL computed using
        %   matrix of predictors X and observed response Y. Data X must be a
        %   numeric matrix of size N-by-P, where P is the number of predictors used
        %   for training this model. Y must be a vector of floating-point numbers
        %   with N elements.
        %
        %   L=LOSS(MODEL,X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
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
        %
        %   See also classreg.learning.regr.RegressionModel,
        %   classreg.learning.regr.RegressionModel/predict.
        
            % Get observation weights
            N = size(X,1);
            args = {                  'lossfun'  'weights'};
            defs = {@classreg.learning.loss.mse  ones(N,1)};
            [funloss,W,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Prepare data
            [X,Y,W] = prepareDataForLoss(this,X,Y,W);
            
            % Check input args
            funloss = classreg.learning.internal.lossCheck(funloss,'regression');
            
            % Get observation margins for hypothesis H.
            Yfit = predict(this,X,extraArgs{:});
            
            % Check
            classreg.learning.internal.regrCheck(Y,Yfit,W);

            % Get loss
            l = funloss(Y,Yfit,W);
        end
     end

end
