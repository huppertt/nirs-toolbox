classdef MajorityUndersampler < classreg.learning.generator.Generator

%   Copyright 2012 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=protected)
        % Logical matrix of size N-by-K for N observations and K classes
        % with class membership.
        C = [];
        
        % Number of observations in the smallest class
        NumSmallest = [];
        
        % Vector of size K for K classes with sampling factors wrt the
        % smallest class
        RatioToSmallest = [];
    end
    
    methods(Hidden)
        function this = MajorityUndersampler(X,Y,W,fitData,classnames,ratioToSmallest)
            this = this@classreg.learning.generator.Generator(X,Y,W,fitData);
            
            this.C = classreg.learning.internal.classCount(classnames,Y);
            K = numel(classnames);
            
            % Find the smallest class with some observations
            sumC = sum(this.C,1);
            sumC(sumC==0) = [];
            if isempty(sumC)
                error(message('stats:classreg:learning:generator:MajorityUndersampler:MajorityUndersampler:AllClassesEmpty'));
            end
            this.NumSmallest = min(sumC);
            
            if isempty(ratioToSmallest) || strcmpi(ratioToSmallest,'default')
                this.RatioToSmallest = ones(1,K);
            else
                if ~isnumeric(ratioToSmallest) || ~isvector(ratioToSmallest) ...
                        || any(ratioToSmallest<0) || all(ratioToSmallest==0)
                    error(message('stats:classreg:learning:generator:MajorityUndersampler:MajorityUndersampler:BadRatioToSmallest'));
                end
                if isscalar(ratioToSmallest)
                    this.RatioToSmallest = repmat(ratioToSmallest,1,K);
                else
                    if numel(ratioToSmallest)~=K ...
                            || any(isnan(ratioToSmallest)) || any(isinf(ratioToSmallest))
                        error(message('stats:classreg:learning:generator:MajorityUndersampler:MajorityUndersampler:RatioToSmallestNaNorInf', K));
                    end
                    this.RatioToSmallest = ratioToSmallest(:)';
                end
            end
        end
    end
    
    methods(Static)
        function [dorus,undersamplerArgs,otherArgs] = processArgs(varargin)
            args = {'ratioToSmallest'};
            defs = {               []};
            [ratioToSmallest,~,otherArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            dorus = ~isempty(ratioToSmallest);
            undersamplerArgs = {};
            if dorus
                undersamplerArgs = {'ratioToSmallest' ratioToSmallest};
            end
        end
    end

    methods
        function [this,X,Y,W,fitData,optArgs] = generate(this)
            % How many samples to make?
            K = size(this.C,2);
            NumPerClass = sum(this.C,1);
            NumToSample = ceil(this.NumSmallest*this.RatioToSmallest);
            NumToSample(NumPerClass==0) = 0;
            
            % Loop over classes and sample the desired number of
            % observations for each class.
            idx = zeros(sum(NumToSample),1);
            idxbegin = 1;
            for k=1:K
                if NumToSample(k)>0
                    idxk = find(this.C(:,k));
                    if NumPerClass(k)~=NumToSample(k)
                        if NumPerClass(k)<NumToSample(k)
                            replaceArgs = {'replace' true};
                        else
                            replaceArgs = {'replace' false};
                        end
                        idxk = datasample(idxk,NumToSample(k),...
                            'weights',this.W(idxk),replaceArgs{:});
                    end
                    idx(idxbegin:idxbegin+NumToSample(k)-1) = idxk;
                    idxbegin = idxbegin + NumToSample(k);
                end
            end
            
            % Shuffle in case the weak learner is sensitive to the order in
            % which observations are found
            idx = idx(randperm(numel(idx)));
            this.LastUseObsForIter = idx;
            
            % Return
            X = this.X(idx,:);
            Y = this.Y(idx);
            W = ones(numel(idx),1); % weights have been used for sampling
            fitData = this.FitData(idx,:);
            optArgs = {};
        end
        
        function this = update(this,X,Y,W,fitData)
            this.X = X;
            this.Y = Y;
            this.W = W;
            this.FitData = fitData;
        end
    end
    
    methods(Static,Hidden)
        function ratioToSmallest = getArgsFromCellstr(varargin)
            args = {'ratioToSmallest'};
            defs = {               []};
            ratioToSmallest = internal.stats.parseArgs(args,defs,varargin{:});
        end
    end
    
end
