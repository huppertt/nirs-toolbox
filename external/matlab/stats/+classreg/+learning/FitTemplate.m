classdef FitTemplate < classreg.learning.internal.DisallowVectorOps

%   Copyright 2010-2014 The MathWorks, Inc.

    
    properties(GetAccess=public,SetAccess=protected,Hidden=true)
        Filled = false; % All parameters filled?
        Method = ''; % Method this template is used for, e.g., Tree
        Type = ''; % classification or regression
        BaseFitObjectArgs = {}; % Cell array of arguments for FullClassificationModel or FullRegressionModel
        MakeModelParams = []; % Function that makes a ModelParams from input args
        MakeFitObject = []; % Function that creates a fit object
        PrepareData = []; % Function to prepare data before fitting
        NamesDataPrepIn = {}; % Names of extra input arguments (besides X,Y) for PrepareData function
        NDataPrepOut = []; % Number of extra output arguments (besides X,Y) from PrepareData function
        MakeModelInputArgs = {}; % Arguments passed to MakeModelParams
        CVPartitionSize = []; % Size of the passed cvpartition object, empty if no object is passed
    end

    properties(GetAccess=public,SetAccess=public,Hidden=true)
        ModelParams = []; % Object of class Params for the fitted model
    end
    
    properties(Constant=true,GetAccess=public,Hidden=true)
        AllowedBaseFitObjectArgs = {'weights' 'predictornames' 'categoricalpredictors' ...
            'responsename' 'responsetransform' 'classnames' 'cost' 'prior' 'scoretransform'};
    end
        
    methods(Static,Hidden)
        function temp = make(method,varargin)
            % Check the type of the required argument
            if ~ischar(method)
                error(message('stats:classreg:learning:FitTemplate:make:BadArgs'));
            end
            
            % Extract type (classification or regression)
            args = {'type'};
            defs = {    ''};
            [usertype,~,modelArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Check usertype
            if ~isempty(usertype)
                usertype = gettype(usertype);
            end

            % Method
            namesclass = classreg.learning.classificationModels();
            namesreg = classreg.learning.regressionModels();
            [tfclass,locclass] = ismember(lower(method),lower(namesclass));
            [tfreg,locreg] = ismember(lower(method),lower(namesreg));
            if ~tfclass && ~tfreg
                error(message('stats:classreg:learning:FitTemplate:make:UnknownMethod', method));
            end
            if     tfclass && tfreg
                method = namesclass{locclass}; % can get it from namesreg too
                type = usertype;
                
                % If type is not passed for an ensemble method, try to
                % figure it out from learner types. This is useful for
                % users who want to type
                %   fitensemble(X,Y,'Subspace',100,'Discriminant')
                % instead of
                %   fitensemble(X,Y,'Subspace',100,'Discriminant','type','classification')
                if isempty(type) && ismember(method,classreg.learning.ensembleModels())
                    [learners,~,~] = internal.stats.parseArgs({'learners'},{},modelArgs{:});
                    if ischar(learners) || isa(learners,'classreg.learning.FitTemplate')
                        learners = {learners};
                    elseif ~iscell(learners)
                        error(message('stats:classreg:learning:FitTemplate:make:BadLearnerTemplates'));
                    end
                    L = numel(learners);
                    % The user can pass several learner templates, and some
                    % of these learners may be appropriate for
                    % classification, some for regression, and some for
                    % both. The ensemble type cannot be determined
                    % unambiguously unless if all learners are appropriate
                    % for one type of learning *only*. For example, in 12a
                    %   t1 = ClassificationDiscriminant.template
                    %   t2 = ClassificationKNN.template
                    %   fitensemble(X,Y,'Subspace',10,{t1 t2})
                    % is going to work because both discriminant and k-NN
                    % can be used for classification only. If you want to
                    % mix discriminant and tree, you have to specify the
                    % ensemble type explicitly:
                    %   t1 = ClassificationDiscriminant.template
                    %   t2 = ClassificationTree.template
                    %   fitensemble(X,Y,'Bag',10,{t1 t2},'type','classification')
                    types = zeros(L,1); % -1 for regression and 1 for classification
                    for l=1:L
                        meth = learners{l};
                        if isa(meth,'classreg.learning.FitTemplate')
                            meth = meth.Method;
                        end
                        isc = ismember(lower(meth),lower(namesclass));
                        isr = ismember(lower(meth),lower(namesreg));
                        if ~isc && ~isr
                            error(message('stats:classreg:learning:FitTemplate:make:UnknownMethod', meth));
                        end
                        types(l) = isc - isr;
                    end
                    if     all(types==1)
                        type = 'classification';
                    elseif all(types==-1)
                        type = 'regression';
                    end
                end
            elseif tfclass
                method = namesclass{locclass};
                type = 'classification';
            else
                method = namesreg{locreg};
                type = 'regression';
            end
            
            % Make sure the type is consistent
            if ~isempty(usertype) && ~strcmp(usertype,type)
                error(message('stats:classreg:learning:FitTemplate:make:UserTypeMismatch', method, usertype));
            end
            
            % Make template
            temp = classreg.learning.FitTemplate(method,modelArgs);
            temp = fillIfNeeded(temp,type);
        end

        function temp = makeFromModelParams(modelParams,varargin)
            if ~isa(modelParams,'classreg.learning.modelparams.ModelParams')
                error(message('stats:classreg:learning:FitTemplate:makeFromModelParams:BadModelParams'));
            end
            method = modelParams.Method;
            type = modelParams.Type;
            args = varargin;
            if isa(modelParams,'classreg.learning.modelparams.EnsembleParams')
                if isa(modelParams.Generator,'classreg.learning.generator.Resampler')
                    args = [args(:)' {'resample' 'on' ...
                        'replace' modelParams.Generator.Replace 'fresample' modelParams.Generator.FResample}];
                elseif isa(modelParams.Generator,'classreg.learning.generator.SubspaceSampler')
                    args = [args(:)' {'npredtosample' modelParams.Generator.NPredToSample}];
                end
                modelParams.Generator = [];
                modelParams.Modifier = [];
                modelParams.Filled = false;
            end
            temp = classreg.learning.FitTemplate(method,args);
            temp.ModelParams = modelParams;
            temp = fillIfNeeded(temp,type);
        end
        
        function catchType(varargin)
            args = {'type'};
            defs = {    []};
            [type,~,~] = internal.stats.parseArgs(args,defs,varargin{:});
            if ~isempty(type)
                error(message('stats:classreg:learning:FitTemplate:catchType:NonEmptyType'));
            end
        end
    end
        
    methods(Access=protected)
        function this = FitTemplate(method,modelArgs)
            this = this@classreg.learning.internal.DisallowVectorOps();
            this.Method = method;
            this.MakeModelInputArgs = modelArgs;
        end
        
        function tf = isfilled(this)
            % Filled already?
            if this.Filled
                tf = true;
                return;
            end
            
            % Check if there are any empty properties besides
            % this.MakeModelInputArgs
            props = {'Type' 'Method' 'BaseFitObjectArgs' 'ModelParams' ...
                'MakeModelParams' 'MakeFitObject' 'PrepareData' ...
                'NamesDataPrepIn' 'NDataPrepOut'};
            tf = false;
            for i=1:length(props)
                if isempty(this.(props{i}))
                    return;
                end
            end
            tf = true;
        end        
    end
    
    methods(Hidden)
        % Make a fitted object using BaseFitObjectArgs stored in this template.
        % The user can override the stored arguments using varargin.
        % Contents of varargin is not stored in this template - you need to
        % override every time you call fit().
        function obj = fit(this,X,Y,varargin)            
            % Make sure the template has been filled
            if isempty(this.Type)
                error(message('stats:classreg:learning:FitTemplate:fit:BadType', this.Method));
            end
            
            % Prepare cell array of arguments from PrepareData
            dataPrepOut = cell(1,this.NDataPrepOut);
            
            % If argument list is empty, use those given to the
            % constructor.
            if isempty(varargin)
                [X,Y,dataPrepOut{1:this.NDataPrepOut}] = ...
                    this.PrepareData(X,Y,this.BaseFitObjectArgs{:});
                    
            % If argument list is not empty, consider this as an attempt to
            % override those supplied to the constructor.
            else
                % Copy default arguments to dataPrepIn
                nDataPrepIn = length(this.NamesDataPrepIn);
                dataPrepIn = cell(1,nDataPrepIn);
                args = this.NamesDataPrepIn;
                defs = repmat({[]},1,nDataPrepIn);
                [dataPrepIn{1:nDataPrepIn}] = ...
                    internal.stats.parseArgs(args,defs,this.BaseFitObjectArgs{:});
                
                % Override dataPrepIn with arguments from varargin
                defs = dataPrepIn;
                [dataPrepIn{1:nDataPrepIn}] = ...
                    internal.stats.parseArgs(args,defs,varargin{:});
                
                % Convert the final dataPrepIn into varargin for
                % this.PrepareData.
                iset = find(~cellfun(@isempty,dataPrepIn));
                nset = length(iset);
                prepArgs = cell(2*nset,1);
                for i=1:nset
                    j = iset(i);
                    prepArgs(2*i-1:2*i) = [this.NamesDataPrepIn(j) dataPrepIn(j)];
                end
                
                % Invoke this.PrepareData
                [X,Y,dataPrepOut{1:this.NDataPrepOut}] = ...
                    this.PrepareData(X,Y,prepArgs{:});
            end
            
            % Check if the data size has been reduced and error if this is not
            % allowed.
            if ~isempty(this.CVPartitionSize) && size(X,1)~=this.CVPartitionSize
                error(message('stats:classreg:learning:FitTemplate:fit:CVPartitionSizeMismatch',...
                    this.CVPartitionSize,size(X,1)));
            end
            
            % Get weights and the rest of the arguments.
            W = dataPrepOut{1};
            fitArgs = dataPrepOut(2:end);
            
            % Fit
            obj = this.MakeFitObject(X,Y,W,this.ModelParams,fitArgs{:});
        end

        function disp(this)
            isLoose = strcmp(get(0,'FormatSpacing'),'loose');
            if strcmp(this.Method,'ByBinaryRegr')
                temp = this.ModelParams.RegressionTemplate;
                if isempty(temp)
                    fprintf('%s\n',getString(message('stats:classreg:learning:FitTemplate:FitTemplateForByBinaryRegr')));
                    if isLoose
                        fprintf('\n');
                    end
                else
                    disp(temp);
                end
                return;
            end
            if isempty(this.Type)
                type = '';
            else
                type = [this.Type ' '];
            end
            fprintf('%s %s%s.\n',getString(message('stats:classreg:learning:FitTemplate:FitTemplateFor')),...
                type,this.Method);
            if     ~isempty(this.ModelParams)
                disp(this.ModelParams);
            elseif ~isempty(this.MakeModelInputArgs)
                N = numel(this.MakeModelInputArgs)/2;
                n=1:N;
                c = this.MakeModelInputArgs(2*n);
                f = this.MakeModelInputArgs(2*n-1);
                s = cell2struct(c(:),f(:),1);
                disp(s);
            else            
                if isLoose
                    fprintf('\n');
                end
            end
        end
        
        function this = fillIfNeeded(this,type)
            % If no type supplied, exit
            if isempty(type)
                return;
            end
            
            % Check type against supplied earlier
            if ~isempty(this.Type)
                if ~strcmp(type,this.Type)
                    error(message('stats:classreg:learning:FitTemplate:fillIfNeeded:IncompatibleType', this.Type));
                end
            end
            this.Type = type;
            
            % Already filled?
            if isfilled(this)
                return;
            end
            
            % Recipe for preparing data
            if strcmp(this.Method,'ByBinaryRegr')
                this.PrepareData = ...
                    @classreg.learning.classif.ClassifByBinaryRegr.prepareData;
                this.NDataPrepOut = 4; % W,dataSummary,classSummary,scoreTransform
                this.NamesDataPrepIn = ...
                    {'weights' 'predictornames' 'categoricalpredictors' 'responsename' 'classnames' 'cost' 'prior' 'scoretransform'};
            elseif strcmp(this.Type,'classification')
                if     strcmp(this.Method,'Tree')
                    this.PrepareData = @ClassificationTree.prepareData;
                elseif strcmp(this.Method,'SVM')
                    this.PrepareData = @ClassificationSVM.prepareData;
                else
                    this.PrepareData = @classreg.learning.classif.FullClassificationModel.prepareData;
                end
                this.NDataPrepOut = 4; % W,dataSummary,classSummary,scoreTransform
                this.NamesDataPrepIn = ...
                    {'weights' 'predictornames' 'categoricalpredictors' 'responsename' 'classnames' 'cost' 'prior' 'scoretransform'};
            else
                this.PrepareData = @classreg.learning.regr.FullRegressionModel.prepareData;
                this.NDataPrepOut = 3; % W,dataSummary,responseTransform
                this.NamesDataPrepIn = {'weights' 'predictornames' 'categoricalpredictors' 'responsename' 'responsetransform'};
            end
            
            % Extract args for cross-validation and resampling.
            [Nfold,partitionArgs,otherArgs,cvpartsize] = ...
                classreg.learning.generator.Partitioner.processArgs(this.MakeModelInputArgs{:});
            docv = ~isempty(Nfold);
            this.CVPartitionSize = cvpartsize;
            [dobag,sampleArgs,otherArgs] = ...
                classreg.learning.generator.Resampler.processArgs(otherArgs{:});
            [dosubspace,subspaceArgs,otherArgs] = ...
                classreg.learning.generator.SubspaceSampler.processArgs(otherArgs{:});
            [dorus,undersamplerArgs,otherArgs] = ...
                classreg.learning.generator.MajorityUndersampler.processArgs(otherArgs{:});
            
            % Are model params supplied as an object?
            args = {'modelparams'};
            defs = {           []};
            [modelParams,~,otherArgs] = ...
                internal.stats.parseArgs(args,defs,otherArgs{:});
            if ~docv && ~isempty(modelParams)
                error(message('stats:classreg:learning:FitTemplate:fillIfNeeded:ModelParamsWithoutCV'));
            end
            
            % Force bagging based on the method name.
            if strcmp(this.Method,'Bag')
                dobag = true;
            end
            
            % Force subspace based on the method name.
            if strcmp(this.Method,'Subspace')
                dosubspace = true;
            end
            
            % Force majority undersampling based on the method name.
            if strcmp(this.Method,'RUSBoost')
                dorus = true;
            end
            
            % Disallow sampling observations and predictors at the same
            % time
            if dobag && dosubspace
                error(message('stats:classreg:learning:FitTemplate:fillIfNeeded:BothBaggingAndSubspaceNotAllowed'));
            end
            
            % Disallow sampling observations and boosting by majority
            % undersampling at the same time
            if dobag && dorus
                error(message('stats:classreg:learning:FitTemplate:fillIfNeeded:BothBaggingAndRUSBoostNotAllowed'));
            end
            
            % Disallow subspace and boosting by majority undersampling at
            % the same time
            if dosubspace && dorus
                error(message('stats:classreg:learning:FitTemplate:fillIfNeeded:BothSubspaceAndRUSBoostNotAllowed'));
            end
            
            % Extract base arguments contained in NamesDataPrepIn
            f = @(x) any(strncmpi(x,this.NamesDataPrepIn,length(x)));
            loc = find(cellfun(f,otherArgs(1:2:end)));
            loc = loc(:)';
            baseArgs = otherArgs(sort([2*loc-1 2*loc]));
            otherArgs([2*loc-1 2*loc]) = [];
            
            % Remove ScoreTransform from learner arguments
            baseArgsForLearner = baseArgs;
            f = @(x) strncmpi(x,'scoretransform',length(x));
            loc = find(cellfun(f,baseArgs(1:2:end)));
            loc = loc(:)';
            baseArgsForLearner([2*loc-1 2*loc]) = []; 

            % Decode method name
            switch this.Method
                case classreg.learning.simpleModels()
                    if     dobag
                        error(message('stats:classreg:learning:FitTemplate:fillIfNeeded:CannotBagSimpleModel', this.Method));
                    elseif dosubspace
                        error(message('stats:classreg:learning:FitTemplate:fillIfNeeded:CannotSubspaceSimpleModel', this.Method));
                    elseif dorus
                        error(message('stats:classreg:learning:FitTemplate:fillIfNeeded:CannotRUSBoostSimpleModel', this.Method));
                    elseif docv
                        this.MakeModelParams = @classreg.learning.modelparams.EnsembleParams.make;
                        if isempty(modelParams)
                            learnerTemplate = ...
                                classreg.learning.FitTemplate.make(this.Method,'type',this.Type,...
                                baseArgsForLearner{:},otherArgs{:});
                        else
                            learnerTemplate = ...
                                classreg.learning.FitTemplate.makeFromModelParams(modelParams,...
                                baseArgsForLearner{:},otherArgs{:});                            
                        end
                        if strcmp(this.Method,'ECOC')
                            meth = 'PartitionedECOC';
                            this.MakeFitObject = ...
                                    @classreg.learning.partition.ClassificationPartitionedECOC;
                        else
                            meth = 'PartitionedModel';
                            if strcmp(this.Type,'classification')
                                this.MakeFitObject = ...
                                    @classreg.learning.partition.ClassificationPartitionedModel;
                            else
                                this.MakeFitObject = ...
                                    @classreg.learning.partition.RegressionPartitionedModel;
                            end
                        end
                        this.MakeModelInputArgs = [baseArgs partitionArgs ...
                            {'method' meth 'learners' learnerTemplate 'nlearn' Nfold}];
                    else
                        this.MakeModelInputArgs = [baseArgs otherArgs];
                        this.MakeModelParams = str2func(['classreg.learning.modelparams.' this.Method 'Params.make']);
                        if strcmp(this.Type,'classification')
                            if     strcmp(this.Method,'ByBinaryRegr')
                                this.MakeFitObject = @classreg.learning.classif.ClassifByBinaryRegr;
                            else
                                this.MakeFitObject = str2func(['Classification' this.Method]);
                            end
                        else
                            this.MakeFitObject = str2func(['Regression' this.Method]);
                        end
                    end
                case classreg.learning.ensembleModels()           
                    this.MakeModelParams = @classreg.learning.modelparams.EnsembleParams.make;                    
                    if     docv
                        f = @(x) strcmpi(x,'nprint');
                        loc = find(cellfun(f,otherArgs(1:2:end)));
                        loc = loc(:)';
                        otherArgs([2*loc-1 2*loc]) = [];
                        if isempty(modelParams)
                            ensembleTemplate = ...
                                classreg.learning.FitTemplate.make(this.Method,'type',this.Type,...
                                subspaceArgs{:},sampleArgs{:},undersamplerArgs{:},...
                                baseArgsForLearner{:},otherArgs{:},'nprint','off');
                        else
                            if ~isempty(sampleArgs)
                                error(message('stats:classreg:learning:FitTemplate:fillIfNeeded:CannotRedefineBagArgs'));
                            end
                            if ~isempty(subspaceArgs)
                                error(message('stats:classreg:learning:FitTemplate:fillIfNeeded:CannotRedefineSubspaceArgs'));
                            end
                            if ~isempty(undersamplerArgs)
                                error(message('stats:classreg:learning:FitTemplate:fillIfNeeded:CannotRedefineRUSBoostArgs'));
                            end
                            ensembleTemplate = ...
                                classreg.learning.FitTemplate.makeFromModelParams(modelParams,...
                                baseArgsForLearner{:},otherArgs{:},'nprint','off');             
                        end                        
                        loc = find(cellfun(f,this.MakeModelInputArgs(1:2:end)),1,'last');
                        if ~isempty(loc)
                            printArgs = [{'nprint'} this.MakeModelInputArgs(2*loc)];
                        else
                            printArgs = {};
                        end
                        this.MakeModelInputArgs = [baseArgs otherArgs partitionArgs printArgs ...
                            {'method' 'PartitionedEnsemble' 'learners' ensembleTemplate ...
                            'printmsg' 'Completed folds: ' 'savetrainable' true 'nlearn' Nfold}];
                        if strcmp(this.Type,'classification')
                            this.MakeFitObject = ...
                                @classreg.learning.partition.ClassificationPartitionedEnsemble;
                        else
                            this.MakeFitObject = ...
                                @classreg.learning.partition.RegressionPartitionedEnsemble;
                        end
                    elseif dobag
                        this.MakeModelInputArgs = [baseArgs otherArgs ...
                            {'method' this.Method} {'resample' 'on'} sampleArgs];
                        if strcmp(this.Type,'classification')
                            this.MakeFitObject = ...
                                @classreg.learning.classif.ClassificationBaggedEnsemble;
                        else
                            this.MakeFitObject = ...
                                @classreg.learning.regr.RegressionBaggedEnsemble;
                        end
                    elseif dosubspace
                        this.MakeModelInputArgs = [baseArgs otherArgs ...
                            {'method' this.Method} {'subspace' true} subspaceArgs];
                        if strcmp(this.Type,'classification')
                            this.MakeFitObject = @classreg.learning.classif.ClassificationEnsemble;
                        else
                            this.MakeFitObject = @classreg.learning.regr.RegressionEnsemble;
                        end
                    elseif dorus
                        this.MakeModelInputArgs = [baseArgs otherArgs ...
                            {'method' this.Method} undersamplerArgs];
                        if strcmp(this.Type,'classification')
                            this.MakeFitObject = @classreg.learning.classif.ClassificationEnsemble;
                        else
                            this.MakeFitObject = @classreg.learning.regr.RegressionEnsemble;
                        end
                    else
                        this.MakeModelInputArgs = [baseArgs otherArgs {'method' this.Method}];
                        if strcmp(this.Type,'classification')
                            this.MakeFitObject = @classreg.learning.classif.ClassificationEnsemble;
                        else
                            this.MakeFitObject = @classreg.learning.regr.RegressionEnsemble;
                        end
                    end
                otherwise
                    error(message('stats:classreg:learning:FitTemplate:fillIfNeeded:UnknownAlgorithm', this.Method));
            end
            
            % Make an object with model parameters
            if isempty(this.ModelParams) % Make from scratch
                [this.ModelParams,baseArgs] = ...
                    this.MakeModelParams(this.Type,this.MakeModelInputArgs{:});
                for n=1:2:numel(baseArgs)
                    if ~ischar(baseArgs{n})
                        error(message('stats:classreg:learning:FitTemplate:fillIfNeeded:BadBaseFitObjectParameterType'));
                    end
                    if ~any(strncmpi(baseArgs{n},this.AllowedBaseFitObjectArgs,length(baseArgs{n})))
                       error(message('stats:classreg:learning:FitTemplate:fillIfNeeded:UnknownBaseFitObjectParameter', baseArgs{ n }));
                    end
                end
                this.BaseFitObjectArgs = baseArgs;
            else % Override the existing content with this.MakeModelInputArgs
                % Do not process valid arguments for cross-validation or
                % resampling at this time. Simply keep them in
                % MakeModelInputArgs.
                %
                % This piece is executed, for example, when the user passes
                % parameters irrelevant to cross-validation to CROSSVAL
                % method:
                %   t = ClassificationTree.fit(X,Y);
                %   crossval(t,'KFold',5,'nprint',1)
                % The 'KFold' parameter is picked before we get to this
                % place, and 'nprint' is passed to setInputArg.
                if isa(this.ModelParams,'classreg.learning.modelparams.EnsembleParams')
                    [~,~,otherArgs] = ...
                        classreg.learning.generator.Resampler.processArgs(this.MakeModelInputArgs{:});
                    [~,~,otherArgs] = ...
                        classreg.learning.generator.Partitioner.processArgs(otherArgs{:});
                    [~,~,otherArgs] = ...
                        classreg.learning.generator.SubspaceSampler.processArgs(otherArgs{:});
                end
                for n=1:(numel(otherArgs)/2)
                    name = otherArgs{2*n-1};
                    value = otherArgs{2*n};
                    if ~strcmpi(name,'method')
                        if ismember(lower(name),this.AllowedBaseFitObjectArgs)
                            this = setBaseArg(this,name,value);
                        else
                            this = setInputArg(this,name,value);
                        end
                    end
                end
            end
            
            % Filled
            this.Filled = true;
        end

        % Set an input argument for this object by setting a property of
        % class ModelParams. This method should be only used by
        % functions that choose optimal values for parameters behind the
        % scene such as, e.g., object factories.
        function this = setInputArg(this,name,value)
            if ~ischar(name)
                error(message('stats:classreg:learning:FitTemplate:setInputArg:ArgNameNotChar'));
            end
            props = properties(this.ModelParams);
            [tf,loc] = ismember(lower(name),lower(props));
            if ~tf
                error(message('stats:classreg:learning:FitTemplate:setInputArg:ArgNameNotFound', name));
            end
            foundprop = props{loc};
            this.ModelParams.(foundprop) = value;
        end
        
        function out = getInputArg(this,name)
            if ~ischar(name)
                error(message('stats:classreg:learning:FitTemplate:getInputArg:ArgNameNotChar'));
            end
            props = properties(this.ModelParams);
            [tf,loc] = ismember(lower(name),lower(props));
            if ~tf
                error(message('stats:classreg:learning:FitTemplate:getInputArg:ArgNameNotFound', name));
            end
            foundprop = props{loc};
            out = this.ModelParams.(foundprop);
        end
        
        function tf = isemptyInputArg(this,name)
            if ~ischar(name)
                error(message('stats:classreg:learning:FitTemplate:isemptyInputArg:ArgNameNotChar'));
            end
            props = properties(this.ModelParams);
            [tf,loc] = ismember(lower(name),lower(props));
            if ~tf
                error(message('stats:classreg:learning:FitTemplate:isemptyInputArg:ArgNameNotFound', name));
            end
            foundprop = props{loc};
            tf = isempty(this.ModelParams.(foundprop));
        end
        
        function this = setBaseArg(this,name,value)
            if ~ischar(name)
                error(message('stats:classreg:learning:FitTemplate:setBaseArg:ArgNameNotChar'));
            end
            f = @(x) strcmpi(name,x);
            loc = find(cellfun(f,this.BaseFitObjectArgs(1:2:end)));
            if isempty(loc) % add this property
                this.BaseFitObjectArgs(end+1:end+2) = {lower(name) value};
            else % override existing property
                loc = 2*loc - 1;
                if loc(1)+1>length(this.BaseFitObjectArgs)
                    error(message('stats:classreg:learning:FitTemplate:setBaseArg:MissingValueInPropertyList', name));
                end
                this.BaseFitObjectArgs{loc(1)+1} = value;
                if ~isscalar(loc)
                    loc = loc(:)';
                    this.BaseFitObjectArgs([loc(2:end) loc(2:end)+1]) = [];
                end
            end
        end
        
        function this = setType(this,type)
            % This has to be matched to what is done in fillIfNeeded()
            this.Type = gettype(type);
            this.Filled = false;
            this.ModelParams = [];
        end
    end
end

function type = gettype(type)
% Type - classification or regression?
if ~ischar(type)
    error(message('stats:classreg:learning:FitTemplate:gettype:BadType'));
end
loc = find(strncmpi(type,{'classification' 'regression'},length(type)));
if isempty(loc)
    error(message('stats:classreg:learning:FitTemplate:gettype:UnknownType'));
end
if loc==1
    type = 'classification';
else
    type = 'regression';
end
end

% Include explicit declarations for Application Compiler. These function
% handles are formed by str2func and cannot be resolved at compilation.

%#function classreg.learning.modelparams.ByBinaryRegrParams.make
%#function classreg.learning.modelparams.DiscriminantParams.make
%#function classreg.learning.modelparams.ECOCParams.make
%#function classreg.learning.modelparams.EnsembleParams.make
%#function classreg.learning.modelparams.KNNParams.make
%#function classreg.learning.modelparams.NaiveBayesParams.make
%#function classreg.learning.modelparams.SVMParams.make
%#function classreg.learning.modelparams.TreeParams.make
