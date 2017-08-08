classdef Predictor < classreg.learning.internal.DisallowVectorOps
%Predictor Super class for all supervised-learning models.

%   Copyright 2010-2014 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=public,Hidden=true)
        Impl = [];
    end
    
    properties(GetAccess=public,SetAccess=protected,Hidden=true)
        DataSummary = struct('PredictorNames',{},'CategoricalPredictors',[],'ResponseName','');
    end

    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        %PREDICTORNAMES Names of predictors used for this model.
        %   The PredictorNames property is a cell array of strings with names of
        %   predictor variables, one name per column of X.
        %
        %   See also classreg.learning.Predictor.
        PredictorNames;
        
        %CATEGORICALPREDICTORS Indices of categorical predictors.
        %   The CategoricalPredictors property is an array with indices of
        %   categorical predictors. The indices are in the range from 1 to the
        %   number of columns in X.        
        %
        %   See also classreg.learning.Predictor.
        CategoricalPredictors;
        
        %RESPONSENAME Name of the response variable.
        %   The ResponseName property is a string with the name of the response
        %   variable Y.
        %
        %   See also classreg.learning.Predictor.
        ResponseName;
    end

    methods(Abstract)
        varargout = predict(this,X,varargin)
    end
    
    methods(Hidden)
        function disp(this)
            internal.stats.displayClassName(this);
            
            % Body display
            s = propsForDisp(this,[]);
            disp(s);
            
            internal.stats.displayMethodsProperties(this);
        end
    end
    
    methods
        function n = get.PredictorNames(this)
            names = this.DataSummary.PredictorNames;
            if isnumeric(names)
                n = classreg.learning.internal.defaultPredictorNames(names);
            else
                n = names;
            end
        end
        
        function c = get.CategoricalPredictors(this)
            c = this.DataSummary.CategoricalPredictors;
        end
        
        function r = get.ResponseName(this)
            r = this.DataSummary.ResponseName;
        end
    end

    methods(Access=protected)
        function this = Predictor(dataSummary)
            this = this@classreg.learning.internal.DisallowVectorOps();
            this.Impl = [];
            this.DataSummary = dataSummary;
        end
        
        function s = propsForDisp(this,s)
            % The 2nd input argument is a struct accumulating fields to be
            % displayed by derived classes.
            if nargin<2 || isempty(s)
                s = struct;
            else
                if ~isstruct(s)
                    error(message('stats:classreg:learning:Predictor:propsForDisp:BadS'));
                end
            end
            s.PredictorNames = this.PredictorNames;
            s.ResponseName = this.ResponseName;
        end
    end

end
