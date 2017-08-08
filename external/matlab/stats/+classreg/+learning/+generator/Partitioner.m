classdef Partitioner < classreg.learning.generator.Generator

%   Copyright 2010-2011 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=protected)
        KFold = 10;
        LastProcessedFold = 0;
        Partition = [];
    end
    
    methods(Hidden)
        function this = Partitioner(X,Y,W,fitData,cvpart,type,kfold,holdout,leaveout)
            this = this@classreg.learning.generator.Generator(X,Y,W,fitData);
            iscvpart = ~isempty(cvpart);
            iskfold = ~isempty(kfold);
            isholdout = ~isempty(holdout);
            isleaveout = (ischar(leaveout) && strcmpi(leaveout,'on')) ...
                || (islogical(leaveout) && leaveout);
            if iscvpart+iskfold+isholdout+isleaveout>1
                error(message('stats:classreg:learning:generator:Partitioner:Partitioner:TooManyCrossvalOptions'));
            end
            if ~iscvpart && ~iskfold && ~isholdout && ~isleaveout
                error(message('stats:classreg:learning:generator:Partitioner:Partitioner:NoCrossvalOptions'));
            end
            N = size(X,1);
            if     iscvpart
                if ~isa(cvpart,'cvpartition')
                    error(message('stats:classreg:learning:generator:Partitioner:Partitioner:BadCVpartition'));
                end
                this.KFold = cvpart.NumTestSets;
                this.Partition = cvpart;
            elseif iskfold
                if ~isnumeric(kfold) || ~isscalar(kfold) || kfold<2
                    error(message('stats:classreg:learning:generator:Partitioner:Partitioner:BadKfold'));
                end
                this.KFold = min(kfold,size(X,1));
                if strcmp(type,'classification')
                    this.Partition = cvpartition(labels(Y),'kfold',kfold);
                else
                    this.Partition = cvpartition(N,'kfold',kfold);
                end
            elseif isholdout
                if ~isnumeric(holdout) || ~isscalar(holdout) || holdout<0 || holdout>1
                    error(message('stats:classreg:learning:generator:Partitioner:Partitioner:BadHoldout'));
                end
                this.KFold = 1;
                if strcmp(type,'classification')
                    this.Partition = cvpartition(labels(Y),'holdout',holdout);
                else
                    this.Partition = cvpartition(N,'holdout',holdout);
                end
            elseif isleaveout
                this.KFold = N;
                this.Partition = cvpartition(N,'leaveout');
            end
        end
    end

    methods(Static)
        function [Nfold,partitionArgs,otherArgs,cvpartsize] = processArgs(varargin)
            % Decode input args
            args = {'cvpartition' 'crossval' 'kfold' 'holdout' 'leaveout'};
            defs = {           []         ''      []        []         ''};
            [cvpart,crossval,kfold,holdout,leaveout,~,otherArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Init
            Nfold = [];
            partitionArgs = {};

            % Process crossval
            if ~isempty(crossval) && ~strcmpi(crossval,'on') && ~strcmpi(crossval,'off')
                error(message('stats:classreg:learning:generator:Partitioner:processArgs:BadCrossval'));
            end
            docv = strcmpi(crossval,'on');
            
            % Process leaveout
            if ~isempty(leaveout) && ~strcmpi(leaveout,'on') && ~strcmpi(leaveout,'off')
                error(message('stats:classreg:learning:generator:Partitioner:processArgs:BadLeaveout'));
            end
            
            % Process args left to right
            iscvpart = ~isempty(cvpart);
            iskfold = ~isempty(kfold);
            isholdout = ~isempty(holdout);
            isleaveout = strcmpi(leaveout,'on');
            cvpartsize = [];
            if iscvpart+iskfold+isholdout+isleaveout>1
                error(message('stats:classreg:learning:generator:Partitioner:processArgs:TooManyCrossvalOptions'));
            end
            if strcmpi(crossval,'off') && iscvpart+iskfold+isholdout+isleaveout>0
                error(message('stats:classreg:learning:generator:Partitioner:processArgs:MismatchCrossvalOpts'));
            end
            if ~iscvpart && ~docv && ~iskfold && ~isholdout && ~isleaveout
                return;
            end
            if iscvpart
                if ~isa(cvpart,'cvpartition')
                    error(message('stats:classreg:learning:generator:Partitioner:processArgs:BadCVpartition'));
                end
                cvpartsize = cvpart.NumObservations;
                Nfold = cvpart.NumTestSets;
                partitionArgs = {'cvpartition' cvpart};
            end
            if ~iscvpart && ~iskfold && ~isholdout && ~isleaveout
                iskfold = true;
                kfold = 10;
            end
            if iskfold
                if ~isnumeric(kfold) || ~isscalar(kfold) || kfold<=1
                    error(message('stats:classreg:learning:generator:Partitioner:processArgs:BadKfold'));
                end
                Nfold = ceil(kfold);
                partitionArgs = {'kfold' Nfold};
            end
            if isholdout
                Nfold = 1;
                partitionArgs = {'holdout' holdout};
            end
            if isleaveout
                if ~strcmpi(leaveout,'on') && ~strcmpi(leaveout,'off')
                    error(message('stats:classreg:learning:generator:Partitioner:processArgs:BadLeaveout'));
                end
                Nfold = 'leaveout';
                partitionArgs = {'leaveout' 'on'};
            end
        end
    end

    methods
        function [this,X,Y,W,fitData,optArgs] = generate(this)
            this.LastProcessedFold = this.LastProcessedFold + 1;
            if this.LastProcessedFold > this.KFold
                this.LastProcessedFold = this.LastProcessedFold - this.KFold;
            end
            idx = training(this.Partition,this.LastProcessedFold);
            X = this.X(idx,:);
            Y = this.Y(idx);
            W = this.W(idx);
            if isempty(this.FitData)
                fitData = [];
            else
                fitData = this.FitData(idx,:);
            end
            this.LastUseObsForIter = idx;
            optArgs = {};
        end
        
        function this = update(this,X,Y,W,fitData)
        end
    end

    methods(Static,Hidden)
        function [cvpart,kfold,holdout,leaveout] = getArgsFromCellstr(varargin)
            args = {'cvpartition' 'kfold' 'holdout' 'leaveout'};
            defs = {           []      []        []      'off'};
            [cvpart,kfold,holdout,leaveout] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
        end
    end
end

