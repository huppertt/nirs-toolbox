classdef SubspaceSampler < classreg.learning.generator.Generator

%   Copyright 2010-2011 The MathWorks, Inc.


    properties(GetAccess=public,SetAccess=protected)
        PredictorNames = {};
        NPredToSample = [];
        Exhaustive = [];
        NumAllCombinations = [];
    end
    
    methods(Hidden)
        function this = SubspaceSampler(X,Y,W,fitData,predictorNames,npredtosample,exhaustive)
            this = this@classreg.learning.generator.Generator(X,Y,W,fitData);
            
            D = size(X,2);
            if isnumeric(predictorNames)
                predictorNames = classreg.learning.internal.defaultPredictorNames(predictorNames);
            end
            if ~iscellstr(predictorNames) || numel(predictorNames)~=D
                error(message('stats:classreg:learning:generator:SubspaceSampler:SubspaceSampler:BadPredictorNames', D));
            end
            this.PredictorNames = predictorNames;

            if ~isnumeric(npredtosample) || ~isscalar(npredtosample) ...
                    || isnan(npredtosample) || npredtosample<=0 || npredtosample>D
                error(message('stats:classreg:learning:generator:SubspaceSampler:SubspaceSampler:BadNPredToSample', D));
            end
            this.NPredToSample = ceil(npredtosample);
            
            if ~islogical(exhaustive) || ~isscalar(exhaustive)
                error(message('stats:classreg:learning:generator:SubspaceSampler:SubspaceSampler:BadExhaustive'));
            end
            this.Exhaustive = exhaustive;
            
            if exhaustive
                this.NumAllCombinations = nchoosek(size(this.X,2),npredtosample);
            end
        end
    end
    
    methods(Static)
        function [dosubspace,subspaceArgs,otherArgs] = processArgs(varargin)
            % Decode input args
            args = {'subspace' 'npredtosample' 'exhaustive'};
            defs = {     false              []        false};
            [dosubspace,npredtosample,exhaustive,~,otherArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Init
            if ~islogical(dosubspace) || ~isscalar(dosubspace)
                error(message('stats:classreg:learning:generator:SubspaceSampler:processArgs:BadSubspace'));
            end
            if ~islogical(exhaustive) || ~isscalar(exhaustive)
                error(message('stats:classreg:learning:generator:SubspaceSampler:processArgs:BadExhaustive'));
            end
            subspaceArgs = {};
            
            % Process args
            if isempty(npredtosample)
                if dosubspace
                    subspaceArgs = {'npredtosample' 1 'exhaustive' exhaustive};
                end
            else
                if ~isnumeric(npredtosample) || ~isscalar(npredtosample) ...
                        || isnan(npredtosample) || npredtosample<=0
                    error(message('stats:classreg:learning:generator:SubspaceSampler:processArgs:BadNPredToSample'));
                end
                dosubspace = true;
                subspaceArgs = {'npredtosample' npredtosample 'exhaustive' exhaustive};
            end            
         end
    end
    
    methods(Access=protected)
        function this = reservePredForIter(this)
            if this.Exhaustive
                D = size(this.X,2);
                nPredToSample = this.NPredToSample;
                nLearn = this.NumAllCombinations;
                if this.MaxT~=nLearn
                    error(message('stats:classreg:learning:generator:SubspaceSampler:reservePredForIter:BadMaxT', nPredToSample, nLearn));
                end
                idx = combnk(1:D,nPredToSample);
                this.PrivUsePredForIter(1:D,1:this.MaxT) = false;
                for n=1:nLearn
                    this.PrivUsePredForIter(idx(n,:),n) = true;
                end
            else
                this = reservePredForIter@classreg.learning.generator.Generator(this);
            end
        end
    end

    methods
        function [this,X,Y,W,fitData,optArgs] = generate(this)
            % Get data size
            [N,D] = size(this.X);
            
            if this.Exhaustive
                % If exhaustive search through all subsets, copy from the
                % prepared UsePredForIter. COMBNK preserves the original
                % ordering of predictors; no need to sort.
                idxpred = find(this.PrivUsePredForIter(:,this.T+1));
                this.LastUsePredForIter = idxpred;
            else
                % Sample features without replacement with updated
                % probabilities
                nUsed = sum(this.UsePredForIter,2);
                nUsed = nUsed - min(nUsed);
                npredtosample = this.NPredToSample;
                beta = npredtosample/D;
                weights = beta.^nUsed;
                                
                % Sample a random subset of predictors without replacement.
                % Use all observations. Make sure predictor indices are
                % sorted in the ascending order to preserve the original
                % order.
                idxpred = sort(datasample(1:D,npredtosample,'replace',false,'weights',weights));

                % Save indices of predictors sampled this time
                this.LastUsePredForIter = idxpred;
            end
            X = this.X(:,idxpred);
            Y = this.Y;
            W = this.W;
            fitData = this.FitData;
            idxobs = 1:N;
            
            % Save indices of observations sampled this time
            this.LastUseObsForIter = idxobs;
            
            % Return predictor names
            optArgs = {'predictornames' this.PredictorNames(idxpred)};
        end
        
        function this = update(this,X,Y,W,fitData)
        end
    end
    
    methods(Static,Hidden)
        function [npredtosample,exhaustive] = getArgsFromCellstr(varargin)
            args = {'npredtosample' 'exhaustive'};
            defs = {[]                     false};
            [npredtosample,exhaustive] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
        end
    end        
end
