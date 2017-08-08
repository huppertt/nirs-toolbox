classdef TreeParams < classreg.learning.modelparams.ModelParams
%TreeParams Decision trees parameters.
%
%   TreeParams properties:
%       AlgCat            - Algorithm to find categorical splits for classification.
%       MaxCat            - Use inexact search if a categorical predictor has more than MaxCat levels.
%       MergeLeaves       - Flag for merging leaves after tree is grown
%       MinParent         - Minimal size of parent node in tree
%       MaxSplits         - Maximal number of splits to impose
%       MinLeaf           - Minimal size of leaf node in tree
%       NSurrogate        - Number of surrogate splits to find
%       NVarToSample      - Number of predictors to select at random for decision split
%       Prune             - Flag for computing the optimal pruning sequence
%       PruneCriterion    - 'error' or 'impurity' for classification and 'mse' for regression
%       QEToler           - Tolerance on mean squared error per tree node (regression only)
%       SplitCriterion    - 'gdi', 'twoing' or 'deviance' for classification and 'mse' for regression
%       Stream            - Random stream for random selection of predictors

%   Copyright 2010-2014 The MathWorks, Inc.

    properties
        SplitCriterion = [];
        MinParent = [];
        MinLeaf = [];
        MaxSplits = [];
        NVarToSample = [];
        MergeLeaves = [];
        Prune = [];
        PruneCriterion = [];
        QEToler = [];
        NSurrogate = [];
        MaxCat = [];
        AlgCat = [];
        Stream = [];
    end

    methods(Access=protected)
        function this = TreeParams(type,splitcrit,minparent,minleaf,...
                maxsplits,nvartosample,mergeleaves,prune,prunecrit,qetoler,...
                nsurrogate,maxcat,algcat,stream)
            this = this@classreg.learning.modelparams.ModelParams('Tree',type);
            this.SplitCriterion = splitcrit;
            this.MinParent = minparent;
            this.MinLeaf = minleaf;
            this.MaxSplits = maxsplits;
            this.NVarToSample = nvartosample;
            this.MergeLeaves = mergeleaves;
            this.Prune = prune;
            this.PruneCriterion = prunecrit;
            this.QEToler = qetoler;
            this.NSurrogate = nsurrogate;
            this.MaxCat = maxcat;
            this.AlgCat = algcat;
            this.Stream = stream;
        end
    end

    methods(Static,Hidden)
        function [holder,extraArgs] = make(type,varargin)
            % Decode input args
            args = {'splitcriterion' ...
                {'minparentsize' 'minparent'}...
                {'minleafsize' 'minleaf'} ...
                'maxnumsplits' ...
                {'numvariablestosample' 'nvartosample'} ...
                'mergeleaves' ...
                'prune' ...
                'prunecriterion' ...
                {'quadraticerrortolerance' 'qetoler'} ...
                'surrogate' ...
                {'maxnumcategories' 'maxcat'} ...
                'algorithmforcategorical' ...
                'stream'};
            defs = {repmat([],1,13)};
            [splitcrit,minparent,minleaf,maxsplits,...
                nvartosample,mergeleaves,prune,prunecrit,qetoler,...
                nsurrogate,maxcat,algcat,stream,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});

            if ~isempty(splitcrit)
                if     ischar(splitcrit)
                    if strcmpi(type,'regression') && ~strcmpi(splitcrit,'mse')
                        error(message('stats:classreg:learning:modelparams:TreeParams:make:BadRegressionSplitCrit'));
                    end
                    if strcmpi(type,'classification') ...
                            && ~any(strncmpi(splitcrit,{'gdi' 'deviance' 'twoing'},length(splitcrit)))
                        error(message('stats:classreg:learning:modelparams:TreeParams:make:BadClassificationSplitCrit'));
                    end
                else
                    if ~iscellstr(splitcrit) || numel(splitcrit)~=2
                        error(message('stats:classreg:learning:modelparams:TreeParams:make:BadSplitCrit'));
                    end
                end
            end
            
            if ~isempty(minparent) && (~isnumeric(minparent) || ~isscalar(minparent) ...
                    || minparent<=0 || isnan(minparent) || isinf(minparent))
                error(message('stats:classreg:learning:modelparams:TreeParams:make:BadMinParent'));
            end
            if ~ischar(minparent)
                minparent = ceil(minparent);
            end
            
            if ~isempty(minleaf) && (~isnumeric(minleaf) || ~isscalar(minleaf) ...
                    || minleaf<=0 || isnan(minleaf) || isinf(minleaf))
                error(message('stats:classreg:learning:modelparams:TreeParams:make:BadMinLeaf'));
            end
            minleaf = ceil(minleaf);

            if ~isempty(maxsplits) && (~isnumeric(maxsplits) || ~isscalar(maxsplits) ...
                    || maxsplits<0 || isnan(maxsplits) || isinf(maxsplits))
                error(message('stats:classreg:learning:modelparams:TreeParams:make:BadMaxSplits'));
            end
            maxsplits = ceil(maxsplits);

            if ~isempty(nvartosample) && ...
                    ~(ischar(nvartosample) && strcmpi(nvartosample,'all')) && ...
                    ~(isnumeric(nvartosample) && isscalar(nvartosample) ...
                    && nvartosample>0 && ~isnan(nvartosample) && ~isinf(nvartosample))
                error(message('stats:classreg:learning:modelparams:TreeParams:make:BadNvartosample'));
            end
            if isnumeric(nvartosample)
                nvartosample = ceil(nvartosample);
            end
            
            if ~isempty(mergeleaves) && (~ischar(mergeleaves) ...
                    || ~ismember(lower(mergeleaves),{'on' 'off'}))
                error(message('stats:classreg:learning:modelparams:TreeParams:make:BadMergeLeaves'));
            end

            if ~isempty(prune) && (~ischar(prune) ...
                    || ~ismember(lower(prune),{'on' 'off'}))
                error(message('stats:classreg:learning:modelparams:TreeParams:make:BadPrune'));
            end
            
            if ~isempty(prunecrit)
                if ischar(prunecrit)
                    if strcmpi(type,'regression') && ~strcmpi(prunecrit,'mse')
                        error(message('stats:classreg:learning:modelparams:TreeParams:make:BadRegressionPruneCrit'));
                    end
                    if strcmpi(type,'classification') ...
                            && ~any(strncmpi(prunecrit,{'error' 'impurity'},length(prunecrit)))
                        error(message('stats:classreg:learning:modelparams:TreeParams:make:BadClassificationPruneCrit'));
                    end
                else
                    error(message('stats:classreg:learning:modelparams:TreeParams:make:BadPruneCrit'));
                end
            end
            
            if ~isempty(qetoler) && (~isfloat(qetoler) || ~isscalar(qetoler) ...
                    || qetoler<=0 || isnan(qetoler) || isinf(qetoler))
                error(message('stats:classreg:learning:modelparams:TreeParams:make:BadQEtoler'));
            end

            if ~isempty(nsurrogate) ...
                    && (~ischar(nsurrogate) || ~ismember(lower(nsurrogate),{'on' 'off' 'all'})) ...
                    && (~isnumeric(nsurrogate) || ~isscalar(nsurrogate) || nsurrogate<0)
                error(message('stats:classreg:learning:modelparams:TreeParams:make:BadSurrogate'));
            end
            if isnumeric(nsurrogate)
                nsurrogate = ceil(nsurrogate);
            end
            
            if ~isempty(maxcat) && (~isnumeric(maxcat) || ~isscalar(maxcat) ...
                    || maxcat<0 || isnan(maxcat) || isinf(maxcat))
                error(message('stats:classreg:learning:modelparams:TreeParams:make:BadMaxcat'));
            end
            maxcat = ceil(maxcat);

            if ~isempty(algcat)
                if strcmpi(type,'regression') 
                    error(message('stats:classreg:learning:modelparams:TreeParams:make:AlgCatForRegression'));
                end
                if ~ischar(algcat)
                    error(message('stats:classreg:learning:modelparams:TreeParams:make:AlgCatNotChar'));
                end
                allowedVals = {'Exact' 'PullLeft' 'PCA' 'OVAbyClass'};
                tf = strncmpi(algcat,allowedVals,length(algcat));
                if sum(tf)~=1
                    error(message('stats:classreg:learning:modelparams:TreeParams:make:AlgCatUnknownValue'));
                end
                algcat = allowedVals{tf};
            end
                        
            % Make argument holder
            holder = classreg.learning.modelparams.TreeParams(type,splitcrit,minparent,...
                minleaf,maxsplits,nvartosample,mergeleaves,prune,prunecrit,qetoler,...
                nsurrogate,maxcat,algcat,stream);
        end
        
        function this = loadobj(obj)
            found = fieldnames(obj);
            
            if ismember('AlgCat',found) && ~isempty(obj.AlgCat)
                algcat = obj.AlgCat;
            else
                algcat = 'auto';
            end
            
            if ismember('MaxCat',found) && ~isempty(obj.MaxCat)
                maxcat = obj.MaxCat;
            else
                maxcat = 10;
            end
            
            if ismember('NSurrogate',found) && ~isempty(obj.NSurrogate)
                nsurrogate = obj.NSurrogate;
            else
                if ismember('Surrogate',found) && strcmpi(obj.Surrogate,'on')
                    nsurrogate = 10;
                else
                    nsurrogate = 0;
                end
            end
            
            stream = [];
            if ismember('Stream',found)
                stream = obj.Stream;
            end

            if ismember('MaxSplits',found) && ~isempty(obj.MaxSplits)
                maxsplits = obj.MaxSplits;
            else
                maxsplits = double(intmax); % cast to size_t in C++
            end
            
            minparent = obj.MinParent;
            if ischar(minparent) && strcmp(minparent,'OneSplit')
                minparent = 10;
                maxsplits = 1;
            end
            
            this = classreg.learning.modelparams.TreeParams(...
                obj.Type,obj.SplitCriterion,minparent,...
                obj.MinLeaf,maxsplits,obj.NVarToSample,obj.MergeLeaves,...
                obj.Prune,obj.PruneCriterion,obj.QEToler,...
                nsurrogate,maxcat,algcat,stream);
        end
    end

    methods(Hidden)
        function this = fillDefaultParams(this,X,Y,W,dataSummary,classSummary) %#ok<INUSD>
            [N,D] = size(X);
            
            if isempty(this.SplitCriterion)
                if strcmpi(this.Type,'classification')
                    this.SplitCriterion = 'gdi';
                else
                    this.SplitCriterion = 'mse';
                end
            end
            
            if isempty(this.MinParent)
                this.MinParent = 10;
            end
            if isempty(this.MinLeaf)
                this.MinLeaf = 1;
            end
            this.MinParent = max(this.MinParent,2*this.MinLeaf);
            
            if isempty(this.MaxSplits) || this.MaxSplits>N-1
                this.MaxSplits = N-1;
            end
            
            if isempty(this.NVarToSample) || ...
                    (~ischar(this.NVarToSample) && this.NVarToSample>=D)
                this.NVarToSample = 'all';
            end
            
            if isempty(this.MergeLeaves)
                this.MergeLeaves = 'on';
            end            
            if isempty(this.Prune)
                this.Prune = 'on';
            end
            
            if isempty(this.PruneCriterion);
                if strcmpi(this.Type,'classification')
                    this.PruneCriterion = 'error';
                else
                    this.PruneCriterion = 'mse';
                end
            end
            
            if strcmpi(this.Type,'regression') && isempty(this.QEToler)
                this.QEToler = 1e-6;
            end
            
            if     isempty(this.NSurrogate)
                this.NSurrogate = 0;
            elseif strcmpi(this.NSurrogate,'on')
                this.NSurrogate = min(D-1,10);
            elseif isnumeric(this.NSurrogate) && this.NSurrogate>D-1
                this.NSurrogate = 'all';
            end
            
            if isempty(this.MaxCat)
                this.MaxCat = 10;
            end
            if isempty(this.AlgCat)
                this.AlgCat = 'auto';
            end
        end
    end

end
