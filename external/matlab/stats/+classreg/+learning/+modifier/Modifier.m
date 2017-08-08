classdef Modifier < classreg.learning.internal.DisallowVectorOps

%   Copyright 2010 The MathWorks, Inc.


    properties(Constant=true,GetAccess=public,Abstract=true)
        FitInfoDescription;
    end
    
    properties(GetAccess=public,SetAccess=protected)
        % Why terminated
        ReasonForTermination = ...
            getString(message('stats:classreg:learning:modifier:Modifier:ReasonForTermination'));
        
        % Number of performed consecutive data modifications
        T = 0; 
        
        % Shrinkage rate
        LearnRate = 1;
    end
    
    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        % Filled fit info
        FitInfo;
    end
    
    properties(GetAccess=protected,SetAccess=protected)
        % Size of FitInfo for one application
        FitInfoSize = [];
        
        % Number of max allowed consecutive data modifications
        MaxT = 0;
        
        % Array with stored info from applications of MODIFY. For example,
        % this could be an array of classification error per hypothesis or
        % an array of classification margin per observation per hypothesis.
        FullFitInfo = [];

        % Terminated early due to fit problems?
        Terminated = false;        
    end

    methods
        function fi = get.FitInfo(this)
            if any(this.FitInfoSize(:))
                fi = reshape(this.FullFitInfo(1:this.T,:),this.T,this.FitInfoSize);
            else
                fi = [];
            end
        end
    end

    methods(Abstract)
        % Modify data (X,Y,W) given hypothesis H (compact classifier)
        % trained on data (X,Y,W) and return the updated object.
        [this,mustTerminate,X,Y,W,fitData] = modify(this,X,Y,W,H,fitData)
        
        % Produce an object for combining predictions from weak
        % hypotheses.
        combiner = makeCombiner(this)
    end
    
    methods(Hidden)
        function [this,mustTerminate,X,Y,W,fitData] = modifyWithT(this,X,Y,W,H,fitData)
            [this,mustTerminate,X,Y,W,fitData] = modify(this,X,Y,W,H,fitData);
            if ~mustTerminate
                this = updateT(this);
            end
        end
        
        function this = reserveFitInfo(this,T)
            T = ceil(T);
            if T<=0
                error(message('stats:classreg:learning:modifier:Modifier:reserveFitInfo:BadT'));
            end
            this.MaxT = this.T + T;
            if any(this.FitInfoSize(:))
                this.FullFitInfo(this.T+1:this.MaxT,:) = zeros(T,this.FitInfoSize);
                this.FullFitInfo(this.MaxT+1:end,:) = [];
            end
        end
    end
        
    methods(Access=protected)
        function this = Modifier(fitInfoSize,learnRate)
            this = this@classreg.learning.internal.DisallowVectorOps();
            if isempty(fitInfoSize) || ~isnumeric(fitInfoSize) ...
                    || ~isvector(fitInfoSize) || any(fitInfoSize<0)
                error(message('stats:classreg:learning:modifier:Modifier:Modifier:BadFitInfoSize'));
            end
            fitInfoSize = ceil(fitInfoSize);
            
            if isempty(learnRate) || ~isnumeric(learnRate) || ~isscalar(learnRate) ...
                    || learnRate<=0 || learnRate>1
                error(message('stats:classreg:learning:modifier:Modifier:Modifier:BadLearnRate'));
            end
            
            this.FitInfoSize = fitInfoSize;
            this.LearnRate = learnRate;
        end
        
        function this = updateT(this)
            this.T = this.T + 1;
            if this.T>this.MaxT
                error(message('stats:classreg:learning:modifier:Modifier:updateT:MaxTExceeded'));
            end
        end
        
        function tf = wasTerminated(this)
            % If terminated abnormally, stop forever
            tf = this.Terminated;
            if tf
                warning(message('stats:classreg:learning:modifier:Modifier:modify:StoppedEarly', class( this ), this.T));
            end
        end
    end

end
