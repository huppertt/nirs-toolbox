classdef EnsembleParams < classreg.learning.modelparams.ModelParams
%EnsembleParams Parameters for ensemble learning.
%
%   EnsembleParams properties:
%       LearnerTemplates      - Templates for learners in ensemble.
%       NLearn                - Number cycles for ensemble learning.
%       LearnRate             - Parameter for incremental shrinkage.
%       MarginPrecision       - Precision of the margin goal for LPBoost or TotalBoost.
%       RobustErrorGoal       - Classification error goal for RobustBoost.
%       RobustMaxMargin       - Maximal allowed margin for RobustBoost.
%       RobustMarginSigma     - Margin spread for RobustBoost.
%       SortLearnersByWeight  - Sort learners by their weights when compacting.
%       NPrint                - Print-out frequency during training.
%       PrintMsg              - Print-out message.
%       Generator             - Data generation for learners in ensemble.
%       Modifier              - Data modification by trained learners in ensemble.
%       SaveTrainable         - Flag for saving trainable learners in ensemble object.
%       DefaultScore          - Default ensemble response for missing classes.
    
%   Copyright 2010-2014 The MathWorks, Inc.

    
    properties
        LearnerTemplates = []; 
        NLearn = [];
        LearnRate = [];
        MarginPrecision = [];
        RobustErrorGoal = [];
        RobustMaxMargin = [];
        RobustMarginSigma = [];
        SortLearnersByWeight = [];
        NPrint = []; 
        PrintMsg = '';
        Generator = [];
        Modifier = [];
        SaveTrainable = [];
        DefaultScore = [];
    end
    
    properties(GetAccess=protected,SetAccess=protected)
        GeneratorArgs = {};
    end
    
    methods(Access=protected)
        function this = EnsembleParams(type,method,learnerTemplates,...
                nlearn,learnRate,marprec,rbeps,rbtheta,rbsigma,...
                sortlearners,nprint,printmsg,saveTrainable,defaultScore,...
                generatorArgs)
            this = this@classreg.learning.modelparams.ModelParams(method,type);
            this.LearnerTemplates = learnerTemplates;
            this.NLearn = nlearn;
            this.LearnRate = learnRate;
            this.MarginPrecision = marprec;
            this.RobustErrorGoal = rbeps;
            this.RobustMaxMargin = rbtheta;
            this.RobustMarginSigma = rbsigma;
            this.SortLearnersByWeight = sortlearners;
            this.NPrint = ceil(nprint);
            this.PrintMsg = printmsg;
            this.SaveTrainable = saveTrainable;
            this.GeneratorArgs = generatorArgs;
            this.DefaultScore = defaultScore;
        end
    end
    
    methods(Static,Hidden)
        function [holder,extraArgs] = make(type,varargin)
            % Decode input args
            args = {'method' 'learners' 'nlearn' 'learnrate' 'marginprecision' ...
                'robusterrorgoal' 'robustmaxmargin' 'robustmarginsigma' ...
                'sortlearners' ...
                'nprint' 'printmsg' 'savetrainable' 'defaultscore'};
            defs = {      ''         {}       []          []                [] ...
                            [] ...
                              []                 []                  [] ...
                      []         ''              []            NaN};
            [method,learnerTemplates,nlearn,learnRate,marprec,...
                robustErrorGoal,robustMaxMargin,robustMarginSigma,...
                sortlearners,nprint,msg,saveTrainable,defaultScore,~,...
                extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            if ~ismember(method,[classreg.learning.ensembleModels() ...
                    {'PartitionedModel' 'PartitionedEnsemble' 'PartitionedECOC'}])
                error(message('stats:classreg:learning:modelparams:EnsembleParams:make:BadMethod', method));
            end
            
            if ~isempty(learnRate) && ...
                    (~isnumeric(learnRate) || ~isscalar(learnRate) ...
                    || learnRate<=0 || learnRate>1)
                error(message('stats:classreg:learning:modelparams:EnsembleParams:make:BadLearnRate'));
            end
            if ~isempty(learnRate) && ...
                    ~ismember(method,{'AdaBoostM1' 'AdaBoostM2' 'AdaBoostMH' ...
                    'LogitBoost' 'GentleBoost' 'LSBoost' 'RUSBoost' ...
                    'PartitionedEnsemble'})
                error(message('stats:classreg:learning:modelparams:EnsembleParams:make:NoLearnRateForAlg'));
            end
            
            if ~isempty(marprec) && (~isnumeric(marprec) || ~isscalar(marprec) ...
                    || marprec<0 || marprec>1)
                error(message('stats:classreg:learning:modelparams:EnsembleParams:make:BadMarginPrecision'));
            end
            if ~isempty(marprec) && ~ismember(method,{'LPBoost' 'TotalBoost' 'PartitionedEnsemble'})
                error(message('stats:classreg:learning:modelparams:EnsembleParams:make:DisallowMarginPrecision', method));
            end

            if ~isempty(robustErrorGoal) && ...
                    (~isnumeric(robustErrorGoal) || ~isscalar(robustErrorGoal) ...
                    || robustErrorGoal<0 || robustErrorGoal>1)
                error(message('stats:classreg:learning:modelparams:EnsembleParams:make:BadRobustErrorGoal'));
            end
            
            if ~isempty(robustMaxMargin) && ...
                    (~isnumeric(robustMaxMargin) || ~isscalar(robustMaxMargin) || robustMaxMargin<0)
                error(message('stats:classreg:learning:modelparams:EnsembleParams:make:BadRobustMaxMargin'));
            end
            
            if ~isempty(robustMarginSigma) && ...
                    (~isnumeric(robustMarginSigma) || ~isscalar(robustMarginSigma) ...
                    || robustMarginSigma<0)
                error(message('stats:classreg:learning:modelparams:EnsembleParams:make:BadRobustMarginSigma'));
            end
            
            if (~isempty(robustErrorGoal) || ~isempty(robustMaxMargin) || ~isempty(robustMarginSigma)) ...
                    && ~ismember(method,{'RobustBoost' 'PartitionedEnsemble'})
                error(message('stats:classreg:learning:modelparams:EnsembleParams:make:BadRobustParams',method));
            end
            
            if ~isempty(sortlearners) 
                if ~strcmpi(sortlearners,'off') && ~strcmpi(sortlearners,'on')
                    error(message('stats:classreg:learning:modelparams:EnsembleParams:make:BadSortLearners'));
                end
                sortlearners = strcmpi(sortlearners,'on');
            end
            
            if ~isempty(nprint) && ~strcmpi(nprint,'off') && ...
                    (~isnumeric(nprint) || ~isscalar(nprint) || nprint<0)
                error(message('stats:classreg:learning:modelparams:EnsembleParams:make:BadNPrint'));
            end
            if isnumeric(nprint)
                nprint = ceil(nprint);
            end
            
            if ~isempty(msg) && ~ischar(msg)
                error(message('stats:classreg:learning:modelparams:EnsembleParams:make:BadMsg'));
            end
            
            if ~isempty(saveTrainable) && ...
                    ~strcmpi(saveTrainable,'on') && ~strcmpi(saveTrainable,'off') ...
                    && (~islogical(saveTrainable) || ~isscalar(saveTrainable))
                error(message('stats:classreg:learning:modelparams:EnsembleParams:make:BadSaveTrainable'));
            end
            
            if ~isempty(defaultScore) && ~isnumeric(defaultScore)
                error(message('stats:classreg:learning:modelparams:EnsembleParams:make:BadDefaultScore'));
            end

            % Store extra args for generator
            % Resampling args are processed only if resample=on. 
            % CV args are processed if at least one of them is filled with
            % non-default values.
            [dobag,sampleArgs,extraArgs] = ...
                classreg.learning.generator.Resampler.processArgs(extraArgs{:});
            [fresample,replace] = ...
                classreg.learning.generator.Resampler.getArgsFromCellstr(sampleArgs{:});
            if dobag
                resample = 'on';
            else
                resample = 'off';
            end
            
            [~,partitionArgs,extraArgs] = ...
                classreg.learning.generator.Partitioner.processArgs(extraArgs{:});
            [cvpart,kfold,holdout,leaveout] = ...
                classreg.learning.generator.Partitioner.getArgsFromCellstr(partitionArgs{:});
            
            [dosubspace,subspaceArgs,extraArgs] = ...
                classreg.learning.generator.SubspaceSampler.processArgs(extraArgs{:});
            [npredtosample,exhaustive] = ...
                classreg.learning.generator.SubspaceSampler.getArgsFromCellstr(subspaceArgs{:});
                        
            if ~ismember(method,{'Subspace' 'PartitionedEnsemble'}) && dosubspace
                error(message('stats:classreg:learning:modelparams:EnsembleParams:make:SubspaceArgsWithoutSubspace'));
            end
            
            [dorus,undersamplerArgs,extraArgs] = ...
                classreg.learning.generator.MajorityUndersampler.processArgs(extraArgs{:});
            ratioToSmallest = ...
                classreg.learning.generator.MajorityUndersampler.getArgsFromCellstr(undersamplerArgs{:});                        

            if ~ismember(method,{'RUSBoost' 'PartitionedEnsemble'}) && dorus
                error(message('stats:classreg:learning:modelparams:EnsembleParams:make:RatioToSmallestWithoutRUSBoost'));
            end
            if strcmp(method,'RUSBoost') && isempty(ratioToSmallest)
                ratioToSmallest = 'default';
            end
            
            generatorArgs = {resample fresample replace ...
                cvpart kfold holdout leaveout ...
                npredtosample exhaustive ratioToSmallest};

            % If one template, allow not to wrap in a cell array.
            if ~iscell(learnerTemplates) ...
                    && ~isa(learnerTemplates,'classreg.learning.FitTemplate') ...
                    && ~ischar(learnerTemplates)
                error(message('stats:classreg:learning:modelparams:EnsembleParams:make:BadLearnerTemplates'));
            end
            if ischar(learnerTemplates)
                learnerTemplates = classreg.learning.FitTemplate.make(learnerTemplates);
            end
            if isa(learnerTemplates,'classreg.learning.FitTemplate')
                learnerTemplates = {learnerTemplates};
            end
            
            % Number of training iterations.
            if ischar(nlearn)
                Ttodo = nlearn;
            else
                Ttodo = nlearn*numel(learnerTemplates);
            end
            
            % Make argument holder
            holder = classreg.learning.modelparams.EnsembleParams(...
                type,method,learnerTemplates,...
                Ttodo,learnRate,marprec,...
                robustErrorGoal,robustMaxMargin,robustMarginSigma,...                
                sortlearners,nprint,msg,saveTrainable,defaultScore,...
                generatorArgs);
        end
    end
    
    
    methods(Access=protected)
        function group = getPropertyGroups(this)
            plist = struct;            

            plist.Type = this.Type;
            plist.Method = this.Method;
            
            str = '';
            for i=1:numel(this.LearnerTemplates)
                temp = this.LearnerTemplates{i};
                if strcmp(temp.Method,'ByBinaryRegr')
                    temp = temp.ModelParams.RegressionTemplate;
                end
                str = [str sprintf('%s',temp.Method)]; %#ok<AGROW>
            end
            plist.LearnerTemplates = str;
            
            plist.NLearn = this.NLearn;
            
            if ismember(this.Method,...
                    {'AdaBoostM1' 'AdaBoostM2' 'AdaBoostMH' ...
                    'LogitBoost' 'GentleBoost' 'LSBoost' 'RUSBoost'})
                plist.LearnRate = this.LearnRate;
            end
            
            if ismember(this.Method,{'LPBoost' 'TotalBoost'})
                plist.MarginPrecision = this.MarginPrecision;
            end            
            
            if strcmp(this.Method,'RobustBoost')
                plist.RobustErrorGoal = this.RobustErrorGoal;
                plist.RobustMaxMargin = this.RobustMaxMargin;
                plist.RobustMarginSigma = this.RobustMarginSigma;
            end
            
            group = matlab.mixin.util.PropertyGroup(plist,'');
        end        
    end
        
    
    methods(Hidden)
        function gen = makeGenerator(this,X,Y,W,fitData,dataSummary,classSummary)
            % Get type (class or reg) and method
            type = this.Type;
            
            % Get generator args
            args = this.GeneratorArgs;
            resample = lower(args{1});
            fresample = args{2};
            replace = lower(args{3});
            cvpart = args{4};
            kfold = args{5};
            holdout = args{6};
            leaveout = lower(args{7});
            
            % Generator arguments above 7 are added in 12a. The checks are
            % for backward compatibility.
            if numel(args)<8
                npredtosample = [];
            else
                npredtosample = args{8};
            end
            if numel(args)<9
                exhaustive = false;
            else
                exhaustive = args{9};
            end
            if numel(args)<10
                ratioToSmallest = [];
            else
                ratioToSmallest = args{10};
            end
            
            %
            % Check generator args
            %
            
            if ~strcmp(resample,'on') && ~strcmp(resample,'off')
                error(message('stats:classreg:learning:modelparams:EnsembleParams:makeGenerator:BadResample'));
            end
            resample = strcmp(resample,'on');
            
            % kfold, holdout and leaveout have been tested by
            % PartitionedModel already
            crossval = ~isempty(cvpart) || ~isempty(kfold) || ...
                ~isempty(holdout) || strcmp(leaveout,'on');
            
            % RUSBoost?
            rusboost = ~isempty(ratioToSmallest);
            
            % Random subspace?            
            subspace = ~isempty(npredtosample);

            % Can do only one of: resampling, cross-validation, subspace
            % sampling or boosting by undersampling. This message verifies
            % the internal consistency of the implementation. Condition
            % resample+crossval+subspace+rusboost>1 should never be
            % satisfied unless there is a coding error in this class or
            % FitTemplate.
            if resample+crossval+subspace+rusboost>1
                error(message('stats:classreg:learning:modelparams:EnsembleParams:makeGenerator:TooManyGeneratorsRequested'));
            end
            
            %
            % Make data generator
            %
            
            if resample
                gen = classreg.learning.generator.Resampler(X,Y,W,fitData,...
                    fresample,replace);
                return;
            end
            
            if crossval
                gen = classreg.learning.generator.Partitioner(X,Y,W,fitData,...
                    cvpart,type,kfold,holdout,leaveout);
                return;
            end
            
            if subspace
                gen = classreg.learning.generator.SubspaceSampler(X,Y,W,fitData,...
                    dataSummary.PredictorNames,npredtosample,exhaustive);
                return;
            end
            
            if rusboost
                gen = classreg.learning.generator.MajorityUndersampler(X,Y,W,...
                    fitData,classSummary.ClassNames,ratioToSmallest);
                return;
            end
            
            % Blank generator by default
            gen = classreg.learning.generator.BlankGenerator(X,Y,W,fitData);
        end
        
        function mod = makeModifier(this,X,Y,W,classSummary)
            if strcmp(this.Type,'classification') ...
                    && any(ismember(this.Method,{'LSBoost'}))
                error(message('stats:classreg:learning:modelparams:EnsembleParams:makeModifier:IncompatibleTypeAndMethod', this.Method, this.Type));
            end
            if strcmp(this.Type,'regression') ...
                    && any(ismember(this.Method,...
                    {'AdaBoostM1' 'AdaBoostM2' 'AdaBoostMH' 'RobustBoost' ...
                    'LogitBoost' 'GentleBoost' 'RUSBoost' 'LPBoost' 'TotalBoost'}))
                error(message('stats:classreg:learning:modelparams:EnsembleParams:makeModifier:IncompatibleTypeAndMethod', this.Method, this.Type));
            end
            switch this.Method
                case 'PartitionedModel'
                    mod = classreg.learning.modifier.BlankModifier();
                case 'PartitionedEnsemble'
                    mod = classreg.learning.modifier.BlankModifier();
                case 'PartitionedECOC'
                    mod = classreg.learning.modifier.BlankModifier();
                case 'AdaBoostM1'
                    mod = classreg.learning.modifier.AdaBoostM1(this.LearnRate);
                case 'AdaBoostM2'
                    mod = classreg.learning.modifier.AdaBoostM2(...
                        classSummary.NonzeroProbClasses,this.LearnRate);
                case 'AdaBoostMH'
                    mod = classreg.learning.modifier.AdaBoostMH(...
                        classSummary.NonzeroProbClasses,this.LearnRate);
                case 'RobustBoost'
                    mod = classreg.learning.modifier.RobustBoost(...
                        this.RobustErrorGoal,this.RobustMaxMargin,this.RobustMarginSigma);
                case 'LogitBoost'
                    mod = classreg.learning.modifier.LogitBoost(this.LearnRate);
                case 'GentleBoost'
                    mod = classreg.learning.modifier.GentleBoost(this.LearnRate);
                case 'LSBoost'
                    mod = classreg.learning.modifier.LSBoost(this.LearnRate);
                case 'RUSBoost'
                    mod = classreg.learning.modifier.RUSBoost(X,Y,W,this.LearnRate);
                case 'LPBoost'
                    mod = classreg.learning.modifier.LPBoost(this.MarginPrecision,numel(W));
                case 'TotalBoost'
                    mod = classreg.learning.modifier.TotalBoost(this.MarginPrecision,numel(W));
                case 'Bag'
                    mod = classreg.learning.modifier.BlankModifier();
                case 'Subspace'
                    mod = classreg.learning.modifier.BlankModifier();
                otherwise
                    error(message('stats:classreg:learning:modelparams:EnsembleParams:makeModifier:UnknownModifier', this.Method));
            end
        end

        function this = fillDefaultParams(this,X,Y,W,dataSummary,classSummary)
            % Check ensemble size
            if     isempty(this.NLearn)
                this.NLearn = numel(this.LearnerTemplates);
            elseif strcmp(this.NLearn,'leaveout')
                % nlearn parameter is set to 'leaveout' by FitTemplate if
                % the user sets 'leaveout' argument to fitensemble to 'on'
                this.NLearn = size(X,1);
            elseif ischar(this.NLearn) && ...
                    strncmpi(this.NLearn,'AllPredictorCombinations',length(this.NLearn)) ...
                    && strcmp(this.Method,'Subspace')
                % nlearn parameter can be set to 'AllPredictorCombinations'
                % by the user for Subspace method
                this.NLearn = 'AllPredictorCombinations'; % work out the numeric value later
                this.GeneratorArgs{9} = true; % set the exhaustive flag
            else
                if ~isnumeric(this.NLearn) || ~isscalar(this.NLearn) || this.NLearn<=0
                    error(message('stats:classreg:learning:modelparams:EnsembleParams:fillDefaultParams:BadNLearn'));
                end
                this.NLearn = ceil(this.NLearn);
            end
            
            % Shrinkage
            if isempty(this.LearnRate)
                this.LearnRate = 1;
            else
                if ~isnumeric(this.LearnRate) || ~isscalar(this.LearnRate) ...
                        || this.LearnRate<=0 || this.LearnRate>1
                    error(message('stats:classreg:learning:modelparams:EnsembleParams:fillDefaultParams:BadLearnRate'));
                end
            end
            
            % Check parameters for LPBoost and TotalBoost
            if ismember(this.Method,{'LPBoost' 'TotalBoost'})
                if isempty(this.MarginPrecision)
                    this.MarginPrecision = 0.01;
                else
                    this.MarginPrecision = max(this.MarginPrecision,1/numel(W));
                end
            end
            
            % Check accuracy parameters for RobustBoost
            if strcmp(this.Method,'RobustBoost')
                if isempty(this.RobustErrorGoal)
                    this.RobustErrorGoal = 0.1;
                end
                if isempty(this.RobustMaxMargin)
                    this.RobustMaxMargin = 0;
                end
                if isempty(this.RobustMarginSigma)
                    this.RobustMarginSigma = 0.1;
                end
            end
            
            % Sort learners by weight?
            if isempty(this.SortLearnersByWeight)
                if ismember(this.Method,{'LPBoost' 'TotalBoost'})
                    this.SortLearnersByWeight = true;
                else
                    this.SortLearnersByWeight = false;
                end
            end
            
            % Check displayed message
            if isempty(this.PrintMsg)
                this.PrintMsg = 'Grown weak learners: ';
            else
                if ~ischar(this.PrintMsg)
                    error(message('stats:classreg:learning:modelparams:EnsembleParams:fillDefaultParams:BadPrintMsg'));
                end
            end
            
            % Save trainable?
            if isempty(this.SaveTrainable)
                this.SaveTrainable = false;
            else
                if ~islogical(this.SaveTrainable)
                    if ~strcmpi(this.SaveTrainable,'on') && ~strcmpi(this.SaveTrainable,'off')
                        error(message('stats:classreg:learning:modelparams:EnsembleParams:fillDefaultParams:BadSaveTrainable'));
                    end
                end
            end

            % Check weak learner templates
            templates = this.LearnerTemplates;
            if isempty(templates)
                error(message('stats:classreg:learning:modelparams:EnsembleParams:fillDefaultParams:NoTemplatesFound'));
            end
            for l=1:length(templates)
                if ~isa(templates{l},'classreg.learning.FitTemplate')
                    error(message('stats:classreg:learning:modelparams:EnsembleParams:fillDefaultParams:NotFitTemplate'));
                end
            end
            if strcmp(this.NLearn,'AllPredictorCombinations') && length(templates)>1
                error(message('stats:classreg:learning:modelparams:EnsembleParams:fillDefaultParams:TooManyLearnersForAllPredictorCombinations'));
            end
            
            % Check classes
            if any(ismember(this.Method,...
                    {'AdaBoostM1' 'RobustBoost' 'LogitBoost' 'GentleBoost'})) ...
                    && length(classSummary.ClassNames)>2
                error(message('stats:classreg:learning:modelparams:EnsembleParams:fillDefaultParams:TooManyClasses', this.Method));
            end
            if any(ismember(this.Method,{'AdaBoostM2'})) ...
                    && length(classSummary.ClassNames)<=2
                error(message('stats:classreg:learning:modelparams:EnsembleParams:fillDefaultParams:TooFewClasses', this.Method));
            end

            % Make special fit info for data generator
            fitData = [];
            switch this.Method
                case {'AdaBoostM2' 'RUSBoost'}
                    % fitData is weights for false hypotheses
                    C = classreg.learning.internal.classCount(...
                        classSummary.NonzeroProbClasses,Y);
                    fitData = repmat(W(:),1,length(classSummary.NonzeroProbClasses));
                    fitData = fitData.*(~C);
                    if any(fitData(:))
                        fitData = fitData / sum(fitData(:));
                    else
                        fitData(:) = 0;
                    end
                case 'AdaBoostMH'
                    % fitData is weights per class
                    fitData = repmat(W(:),1,length(classSummary.NonzeroProbClasses));
                    if any(fitData(:))
                        fitData = fitData / sum(fitData(:));
                    else
                        fitData(:) = 0;
                    end
                case 'RobustBoost'
                    % fitData is current margins
                    fitData = zeros(numel(Y),1);
                case 'LogitBoost'
                    % fitData has 2 columns.
                    % 1st column is 0/1 class label.
                    % 2nd column is current predicted score for the 1st
                    % class in the range from -Inf to +Inf.
                    fitData = zeros(numel(Y),2);
                    fitData(:,1) = classSummary.NonzeroProbClasses(1)==Y;
                    Y = 4*fitData(:,1)-2;
                case 'GentleBoost'
                    % fitData is a column-vector with accumulated model
                    % predictions for training data
                    fitData = zeros(numel(Y),1);
                    Y = double(classSummary.NonzeroProbClasses(1)==Y);
                    Y(Y==0) = -1;
                case 'TotalBoost'
                    % fitData are the original observation weights
                    fitData = W(:);
            end
            
            % Set default score for missing classes
            if ismember(this.Method,...
                    {'AdaBoostM1' 'AdaBoostM2' 'AdaBoostMH' 'RobustBoost' 'LSBoost' ...
                    'GentleBoost' 'LogitBoost' 'RUSBoost' 'LPBoost' 'TotalBoost'})
                this.DefaultScore = -Inf;
            end
            
            % Make data generator
            if isempty(this.Generator)
                this.Generator = makeGenerator(this,X,Y,W,fitData,dataSummary,classSummary);
            end
            
            % After the generator is made, get the number of combinations
            % for exhaustive search through all subspaces of fixed
            % dimensionality
            if strcmp(this.NLearn,'AllPredictorCombinations')
                this.NLearn = this.Generator.NumAllCombinations;
            end
            
            % Make data modifier
            if isempty(this.Modifier)
                this.Modifier = makeModifier(this,X,Y,W,classSummary);
            end

            % Loop through weak learner templates and reset their params if
            % needed
            for l=1:length(templates)
                learner = templates{l};
                
                % Special case: The learner template needs to be wrapped in
                % a ByBinaryRegr class.
                if any(ismember(this.Method,{'LogitBoost' 'GentleBoost'}))
                    if ~ismember(learner.Method,[{'ByBinaryRegr'} classreg.learning.regressionModels()])
                        error(message('stats:classreg:learning:modelparams:EnsembleParams:fillDefaultParams:NoRegressionLearner',learner.Method));
                    end
                    if ~strcmp(learner.Method,'ByBinaryRegr')
                        learner = setType(learner,'regression');
                        learner = fillIfNeeded(learner,'regression');
                        learner = setBaseArg(learner,'predictornames',dataSummary.PredictorNames);
                        learner = setBaseArg(learner,'responsename',dataSummary.ResponseName);
                        learner = setBaseArg(learner,'categoricalpredictors',dataSummary.CategoricalPredictors);
                        learner = classreg.learning.FitTemplate.make('ByBinaryRegr','learner',learner);
                        learner = setBaseArg(learner,'classnames',classSummary.NonzeroProbClasses);
                        learner = setBaseArg(learner,'prior','empirical');
                        K = length(classSummary.NonzeroProbClasses);
                        learner = setBaseArg(learner,'cost',ones(K)-eye(K));
                    end
                end
                
                % Fill learner template for correct type (classification or
                % regression)
                learner = fillIfNeeded(learner,this.Type);
                
                % Copy class summary from the ensemble or partitioned model
                % into the weak learner.
                if strcmp(this.Type,'classification') && ~strcmp(learner.Method,'ByBinaryRegr')
                    if strcmp(this.Method,'PartitionedECOC') ...
                            || strcmp(learner.Method,'NaiveBayes')
                        % Pass all classes to each fold for ECOC and
                        % cross-validated NaiveBayes models including
                        % classes with zero probabilities. For ECOC, the
                        % user can pass a custom coding matrix. For
                        % NaiveBayes, the user can specify a matrix of
                        % kernel widths with one row per class. In either
                        % case, the number and relative class positions in
                        % the ClassNames property cannot change across
                        % folds.
                        K = length(classSummary.ClassNames);
                        learner = setBaseArg(learner,'classnames',classSummary.ClassNames);
                        [~,pos] = ismember(classSummary.NonzeroProbClasses,classSummary.ClassNames);
                        prior = zeros(1,K);
                        prior(pos) = classSummary.Prior;
                        learner = setBaseArg(learner,'prior',prior);
                        if ~isempty(classSummary.Cost)
                            cost = zeros(K);
                            cost(pos,pos) = classSummary.Cost;
                            learner = setBaseArg(learner,'cost',cost);
                        end
                    elseif ismember(this.Method,{'PartitionedModel' 'PartitionedEnsemble'})
                        % Enforce prior and cost specified by the user or
                        % assumed by default for the entire data in each
                        % partition. This corresponds to cross-validation
                        % with stratification by class. Use only classes
                        % with non-zero probabilities in each fold.
                        learner = setBaseArg(learner,'classnames',classSummary.NonzeroProbClasses);
                        learner = setBaseArg(learner,'prior',classSummary.Prior);
                        learner = setBaseArg(learner,'cost',classSummary.Cost);
                    else
                        % The prior is always set to 'empirical' - this is
                        % just an input argument to prepareData() to be run
                        % for every weak learner. The prior supplied to the
                        % ensemble is used to set the initial observation
                        % weights. For boosting, setting prior for weak
                        % learners to the ensemble prior would screw up
                        % weight updates. For bagging, this would change
                        % bootstrap replicas.
                        %
                        % The cost matrix should be always set to the
                        % default. Costs for ensembles are included by
                        % adjusting prior probabilities.
                        learner = setBaseArg(learner,'classnames',classSummary.NonzeroProbClasses);
                        learner = setBaseArg(learner,'prior','empirical');
                        K = length(classSummary.NonzeroProbClasses);
                        learner = setBaseArg(learner,'cost',ones(K)-eye(K));
                    end
                end
                
                % Copy data summary from the ensemble into the weak
                % learner. Assume that Generator is not allowed to
                % sub-sample input predictors and so the dimensionality of
                % the returned data is always the same. Therefore the list
                % of predictor names and indices of categorical predictors
                % are fixed.
                learner = setBaseArg(learner,'predictornames',dataSummary.PredictorNames);
                learner = setBaseArg(learner,'categoricalpredictors',dataSummary.CategoricalPredictors);
                learner = setBaseArg(learner,'responsename',dataSummary.ResponseName);

                % Set up arguments for specific learners
                switch learner.Method
                    case 'Tree'
                        if strcmp(this.Method,'PartitionedModel')
                            this.DefaultScore = 0;
                        end
                        if ~strcmp(this.Method,'PartitionedModel')
                            if isempty(learner.ModelParams.MergeLeaves)
                                learner.ModelParams.MergeLeaves = 'off';
                            end
                            if isempty(learner.ModelParams.Prune)
                                learner.ModelParams.Prune = 'off';
                            end
                        end
                        if ismember(this.Method,{'AdaBoostM1' 'AdaBoostM2' ...
                                'AdaBoostMH' 'RobustBoost' ...
                                'RUSBoost' 'LPBoost' 'TotalBoost'})
                            if isempty(learner.ModelParams.MinParent) ...
                                    && isempty(learner.ModelParams.MinLeaf)
                                learner.ModelParams.MinLeaf   = 1;
                                learner.ModelParams.MinParent = 2;
                                if isempty(learner.ModelParams.MaxSplits)
                                    learner.ModelParams.MaxSplits = 1;
                                end
                            end
                        end
                        if strcmp(this.Method,'LSBoost')
                            if isempty(learner.ModelParams.MinParent) ...
                                    && isempty(learner.ModelParams.MinLeaf)
                                learner.ModelParams.MinLeaf   = 5;
                                learner.ModelParams.MinParent = 10;
                                if isempty(learner.ModelParams.MaxSplits)
                                    learner.ModelParams.MaxSplits = 1;
                                end
                            end
                        end
                        if strcmp(this.Method,'Bag')
                            this.DefaultScore = 0;
                            if isempty(learner.ModelParams.NVarToSample)
                                p = size(X,2);
                                if strcmp(this.Type,'classification')
                                    learner.ModelParams.NVarToSample = ceil(sqrt(p));
                                else
                                    learner.ModelParams.NVarToSample = ceil(p/3);
                                end
                            end
                            if isempty(learner.ModelParams.MinParent) ...
                                    && isempty(learner.ModelParams.MinLeaf)
                                if strcmp(this.Type,'classification')
                                    learner.ModelParams.MinLeaf   = 1;
                                    learner.ModelParams.MinParent = 2;
                                else
                                    learner.ModelParams.MinLeaf   =  5;
                                    learner.ModelParams.MinParent = 10;
                                end
                            end
                        end
                        if isa(this.Generator,'classreg.learning.generator.SubspaceSampler')
                            error(message('stats:classreg:learning:modelparams:EnsembleParams:fillDefaultParams:TreesNotAllowedForSubspace'));
                        end
                        switch this.Method
                            case {'AdaBoostM1' 'AdaBoostMH' 'RobustBoost'}
                                learner = setBaseArg(learner,...
                                    'ScoreTransform',@classreg.learning.transform.symmetricismax);
                            case {'LPBoost' 'TotalBoost'}
                                learner = setBaseArg(learner,...
                                    'ScoreTransform',@classreg.learning.transform.symmetric);
                        end
                        %{
                    case 'NeuralNet'
                        switch this.Method
                            case 'AdaBoostM1'
                                learner = ...
                                    setBaseArg(learner,'ScoreTransform',@classreg.learning.transform.symmetricismax);
                        end
                        learner = setInputArg(learner,'ShowWindow',false);
                        if isemptyInputArg(learner,'NEpochs')
                            learner = setInputArg(learner,'NEpochs',100);
                        end
                        %}
                    case 'ByBinaryRegr'
                        if strcmp(learner.ModelParams.RegressionTemplate.Method,'Tree')
                            if isempty(learner.ModelParams.RegressionTemplate.ModelParams.MergeLeaves)
                                learner.ModelParams.RegressionTemplate.ModelParams.MergeLeaves = 'off';
                            end
                            if isempty(learner.ModelParams.RegressionTemplate.ModelParams.Prune)
                                learner.ModelParams.RegressionTemplate.ModelParams.Prune = 'off';
                            end
                            if any(ismember(this.Method,{'LogitBoost' 'GentleBoost'}))
                                if isempty(learner.ModelParams.RegressionTemplate.ModelParams.MinParent) ...
                                        && isempty(learner.ModelParams.RegressionTemplate.ModelParams.MinLeaf) ...
                                        && isempty(learner.ModelParams.RegressionTemplate.ModelParams.MaxSplits)
                                    learner.ModelParams.RegressionTemplate.ModelParams.MaxSplits = 1;
                                end
                            end
                        end
                        %{
                        if strcmp(learner.ModelParams.RegressionTemplate.Method,'NeuralNet') ...
                                && any(ismember(this.Method,ensembleModels()))
                            if isemptyInputArg(learner.ModelParams.RegressionTemplate,'ShowWindow')
                                learner.ModelParams.RegressionTemplate = ...
                                    setInputArg(learner.ModelParams.RegressionTemplate,...
                                    'ShowWindow',false);
                            end
                            if isemptyInputArg(learner.ModelParams.RegressionTemplate,'NEpochs')
                                learner.ModelParams.RegressionTemplate = ...
                                    setInputArg(learner.ModelParams.RegressionTemplate,...
                                    'NEpochs',100);
                            end
                        end
                        %}
                    case 'Discriminant'
                        if isempty(learner.ModelParams.FillCoeffs)
                            learner.ModelParams.FillCoeffs = false;
                        end
                        if any(ismember(this.Method,{'PartitionedModel' 'Bag'}))
                            this.DefaultScore = 0; % DA scores are posterior probs
                        end
                        switch this.Method
                            case {'AdaBoostM1' 'AdaBoostMH' 'RobustBoost'}
                                learner = setBaseArg(learner,...
                                    'ScoreTransform',@classreg.learning.transform.symmetricismax);
                            case {'LPBoost' 'TotalBoost'}
                                learner = setBaseArg(learner,...
                                    'ScoreTransform',@classreg.learning.transform.symmetric);
                        end
                        if any(ismember(this.Method,classreg.learning.ensembleModels()))
                            if isempty(learner.ModelParams.DiscrimType)
                                learner.ModelParams.DiscrimType = 'pseudoLinear';
                            end                            
                        end
                    case 'KNN'
                        if strcmp(this.Method,'PartitionedModel')
                            this.DefaultScore = 0;
                        elseif ~strcmp(this.Method,'Subspace')
                            error(message('stats:classreg:learning:modelparams:EnsembleParams:fillDefaultParams:UseKNNforSubspaceOnly'));
                        end
                    case 'NaiveBayes'
                        if strcmp(this.Method,'PartitionedModel')
                            this.DefaultScore = 0;
                        else
                            error(message('stats:classreg:learning:modelparams:EnsembleParams:fillDefaultParams:MethodNotAllowedForEnsembleLearning',...
                                learner.Method));
                        end
                    case 'SVM'
                        if strcmp(this.Method,'PartitionedModel')
                            if ~isempty(learner.ModelParams.Alpha)
                                error(message('stats:classreg:learning:modelparams:EnsembleParams:fillDefaultParams:NoSVMAlphaForCrossValidation'));
                            end                            
                            this.DefaultScore = 0;
                        else
                            error(message('stats:classreg:learning:modelparams:EnsembleParams:fillDefaultParams:MethodNotAllowedForEnsembleLearning',...
                                learner.Method));
                        end
%                         switch this.Method
%                             case {'AdaBoostM1' 'AdaBoostMH' 'RobustBoost'}
%                                 learner = setBaseArg(learner,...
%                                     'ScoreTransform',@classreg.learning.transform.symmetricismax);
%                             case {'LPBoost' 'TotalBoost'}
%                                 learner = setBaseArg(learner,...
%                                     'ScoreTransform',@classreg.learning.transform.symmetriclogit);
%                         end
                    case classreg.learning.ensembleModels()
                        if ~strcmp(this.Method,'PartitionedEnsemble')
                            error(message('stats:classreg:learning:modelparams:EnsembleParams:fillDefaultParams:NoEnsembleOfEnsembles'));
                        end
                        if ismember(learner.Method,{'Bag'})
                            this.DefaultScore = 0;
                        else
                            this.DefaultScore = -Inf;
                        end
                    case classreg.learning.weakLearners()
                    %
                    % Append more models here.
                    %
                    otherwise
                        error(message('stats:classreg:learning:modelparams:EnsembleParams:fillDefaultParams:UnknownLearnerMethod', learner.Method));
                end
                
                % Copy updated learner back into the template list
                templates{l} = learner;
            end
            
            % Save the modified weak learner templates
            this.LearnerTemplates = templates;
        end
    end
    
end

