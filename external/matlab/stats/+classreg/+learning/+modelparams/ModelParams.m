classdef ModelParams < classreg.learning.internal.DisallowVectorOps ...
        & matlab.mixin.CustomDisplay
%ModelParams Super class for learning model parameters (before training).
    
%   Copyright 2010-2013 The MathWorks, Inc.

    properties(SetAccess=protected,GetAccess=public)
        Method = ''; % name of the method such as, for example, 'Tree'
        Type = ''; % classification or regression
    end
    
    properties(SetAccess=public,GetAccess=public,Hidden=true)
        Filled = false;% Have all arguments been filled?
    end

    methods(Abstract,Static,Hidden)
        [holder,extraArgs] = make(type,varargin)
    end

    methods(Abstract,Hidden)
        this = fillDefaultParams(this,X,Y,W,dataSummary,classSummary)
    end

    methods(Access=protected)
        function header = getHeader(this)
            header = '';
        end
        
        function this = ModelParams(method,type)
            this = this@classreg.learning.internal.DisallowVectorOps();
            this.Method = method;
            this.Type = type;
        end
    end

    methods(Hidden)
%         function this = fillIfNeeded(this,X,Y,W,dataSummary,classSummary)
%             % original
%             if ~isfilled(this)
%                 this = fillDefaultParams(this,X,Y,W,dataSummary,classSummary);
%             end
%             this.Filled = true;
%         end

        function this = fillIfNeeded(this,X,Y,W,dataSummary,classSummary)
            % Modified 2/28/14
            if ~this.Filled
                this = fillDefaultParams(this,X,Y,W,dataSummary,classSummary);
            end
            this.Filled = true;
        end
        
        function tf = isfilled(this)
            % Filled already?
            if this.Filled
                tf = true;
                return;
            end
            
            % Check if there are any empty properties
            props = properties(this);
            tf = false;
            for i=1:length(props)
                if isempty(this.(props{i}))
                    return;
                end
            end
            tf = true;
        end
    end

end
