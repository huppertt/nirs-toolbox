classdef Generator < classreg.learning.internal.DisallowVectorOps

%   Copyright 2010-2011 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=protected)
        X = [];
        Y = [];
        W = [];
        FitData = [];
        
        % Number of performed consecutive data generations
        T = 0;
    end
    
    properties(GetAccess=public,SetAccess=protected,Dependent=true)
        % Logical N-by-T: true if observation n (n=1,...,N) 
        % is included in iteration t (t=1,...,T) and false otherwise
        UseObsForIter;
        
        % Logical D-by-T: true if predictor d (d=1,...,D) is included in
        % iteration t (t=1,...,T) and false otherwise
        UsePredForIter;
    end
    
    properties(GetAccess=protected,SetAccess=protected)        
        % Number of max allowed consecutive data generations
        MaxT = 0;
        
        % Logical array. True if this observation was used at this iteration.
        PrivUseObsForIter = false(0);
        
        % Indices of observations used at the last iteration.
        LastUseObsForIter = [];

        % Logical array. True if this observation was used at this iteration.
        PrivUsePredForIter = false(0);
        
        % Indices of observations used at the last iteration.
        LastUsePredForIter = [];
    end

    methods(Access=protected)
        function this = Generator(X,Y,W,fitData)
            this = this@classreg.learning.internal.DisallowVectorOps();
            if ~isnumeric(X) || ndims(X)~=2
                error(message('stats:classreg:learning:generator:Generator:Generator:BadX'));
            end
            [N,D] = size(X);
            if (~isnumeric(Y) && ~isa(Y,'classreg.learning.internal.ClassLabel')) ...
                    || ~isvector(Y) || numel(Y)~=N
                error(message('stats:classreg:learning:generator:Generator:Generator:BadY', N));
            end
            if ~isvector(W) || numel(W)~=N
                error(message('stats:classreg:learning:generator:Generator:Generator:BadW', N));
            end
            if ~isempty(fitData) && (~isnumeric(fitData) || size(fitData,1)~=N)
                error(message('stats:classreg:learning:generator:Generator:Generator:BadFitData', N));
            end
            this.X = X;
            this.Y = Y;
            this.W = W;
            this.FitData = fitData;
            this.LastUseObsForIter = 1:N;
            this.LastUsePredForIter = 1:D;
        end
        
        function this = updateT(this)
            this.T = this.T + 1;
            if this.T>this.MaxT
                error(message('stats:classreg:learning:generator:Generator:updateT:MaxTExceeded'));                
            end
            this.PrivUseObsForIter(this.LastUseObsForIter,this.T) = true;            
            this.PrivUsePredForIter(this.LastUsePredForIter,this.T) = true;            
        end
        
        function this = reservePredForIter(this)
            D = size(this.X,2);
            this.PrivUsePredForIter(1:D,this.T+1:this.MaxT) = false;
            this.PrivUsePredForIter(:,this.MaxT+1:end) = [];
        end
    end

    methods(Abstract)
        % Generate data and record what data are generated.
        % X, Y, W, and fitData are generated from the base class
        % properties. 
        % optArgs is a cellstr of optional arguments to be passed to
        % ensemble learners.
        [this,X,Y,W,fitData,optArgs] = generate(this)
        
        % Update the state of the generator given new data.
        this = update(this,X,Y,W,fitData)
    end
    
    methods
        function usenfort = get.UseObsForIter(this)
            usenfort = this.PrivUseObsForIter(:,1:this.T);
        end

        function usenfort = get.UsePredForIter(this)
            if isempty(this.PrivUsePredForIter)
                % PrivUsePredForIter is added in 12a. This branch is for
                % backward compatibility.
                usenfort = true(size(this.X,2),this.T);
            else
                usenfort = this.PrivUsePredForIter(:,1:this.T);
            end
        end
    end
    
    methods(Hidden)
        function this = updateWithT(this,X,Y,W,fitData)
            this = update(this,X,Y,W,fitData);
            this = updateT(this);
        end
        
        function this = reserveFitInfo(this,T)
            N = size(this.X,1);
            T = ceil(T);
            if T<=0
                error(message('stats:classreg:learning:generator:Generator:reserveFitInfo:BadT'));
            end
            this.MaxT = this.T + T;
            this.PrivUseObsForIter(1:N,this.T+1:this.MaxT) = false;
            this.PrivUseObsForIter(:,this.MaxT+1:end) = [];
            this = reservePredForIter(this);
        end
    end

end

