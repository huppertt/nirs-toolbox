classdef NaiveBayesParams < classreg.learning.modelparams.ModelParams
%NaiveBayesParams Parameters for naive Bayes Classification models.
%
%   NaiveBayesParams properties:
%       DistributionNames       - Distribution names.
%       Kernel                  - The type of kernel smoother used for each predictor using the 'kernel' distribution.
%       Support                 - The kernel density support regions. 
%       Width                   - The width of the kernel smoothing window.  

%   Copyright 2013-2014 The MathWorks, Inc.

    properties
        DistributionNames = []; 
        Kernel = []; 
        Support = [];   
        Width = [];       
    end
    
    properties (Constant, Access=protected)
        LegalDistNames = {'kernel','mvmn','normal'};        % Univariate distributions only, so 'mn' special case omitted.
        LegalSupportNames = {'unbounded','positive'};
        LegalKernelNames = {'box','epanechnikov','normal','triangle'};
        LegalCatDistNames = {'mvmn'};
        DefaultCatDistName = 'mvmn';
        DefaultNumDistName = 'normal';
        DefaultKernel = 'normal';
        DefaultSupport = 'unbounded';
        DefaultWidth = NaN;                                 % NaNs are replaced on the fly during fitting.
    end
    
    methods (Access = protected)
        function this = NaiveBayesParams (DistributionNames, Kernel, Support, Width)
            this = this@classreg.learning.modelparams.ModelParams('NaiveBayes','classification');
            this.DistributionNames = DistributionNames(:)'; % Make sure it's a row.
            this.Kernel = Kernel;
            this.Support = Support;
            this.Width = Width;
        end
    end
    
    methods (Static,Hidden)
        function [holder,extraArgs] = make (~,varargin)
            % Decode input args
            args = {'DistributionNames', 'Kernel', 'Support', 'Width'};
            defs =  cell(1,length(args));
            [distributionNames, kernel, support, width, ~, extraArgs] ...
                = internal.stats.parseArgs(args, defs, varargin{:});
            
            % Make argument holder
            holder = classreg.learning.modelparams.NaiveBayesParams(...
                distributionNames, kernel, support, width);
        end
    end
    
    methods (Hidden)
        function this = fillDefaultParams (this,~,~,~,dataSummary,classSummary)
            this = checkDistributionNames(this, numDims(dataSummary), dataSummary);
            this = checkKernelArgs(this, numClasses(classSummary), numDims(dataSummary));
        end
    end
    
    methods (Access=protected)
        function this = checkDistributionNames (this, NumDims, dataSummary)
            if isempty(this.DistributionNames)
                if any(dataSummary.CategoricalPredictors)
                    this.DistributionNames = repmat({this.DefaultNumDistName},1,NumDims);
                    this.DistributionNames(dataSummary.CategoricalPredictors) = {this.DefaultCatDistName};
                else
                    this.DistributionNames = this.DefaultNumDistName;
                end
            end
            if ischar(this.DistributionNames)
                if strcmpi(this.DistributionNames, 'mn')
                    if any(dataSummary.CategoricalPredictors)
                        error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:MnAndCategorical'));
                    end
                    this.DistributionNames = 'mn';
                else
                    [this.DistributionNames, s] = checkAndCompleteString(this.DistributionNames, this.LegalDistNames, 3);
                    if isempty(this.DistributionNames)
                        error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:DistributionNameUnknown', s));
                    end
                end
                % Verify that CategoricalPredictors have legal distributions
                DistCell = repmat({this.DistributionNames}, 1, NumDims);
                if ~all(strcmpi(DistCell(dataSummary.CategoricalPredictors), this.LegalCatDistNames))
                    error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:DistributionNameNotCat'));
                end
            elseif iscell(this.DistributionNames)
                if length(this.DistributionNames) ~= NumDims
                    error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:DistributionNameBadLength'));
                end
                if ~all(cellfun(@ischar, this.DistributionNames))
                    error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:DistributionNameBadType'));
                end
                if any(strcmpi('mn', this.DistributionNames))
                    error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:MnInCellarray'));
                end
                for d=1:NumDims
                    [this.DistributionNames{d}, s] = checkAndCompleteString(this.DistributionNames{d}, this.LegalDistNames, 3);
                    if isempty(this.DistributionNames{d})
                        error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:DistributionNameUnknown', s));
                    end
                end
                % Verify that CategoricalPredictors have legal distributions
                if ~all(strcmpi(this.DistributionNames(dataSummary.CategoricalPredictors), this.LegalCatDistNames))
                    error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:DistributionNameNotCat'));
                end
            else
                error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:DistributionNameBadType'));
            end
        end
        
        function this = checkKernelArgs (this, NumClasses, NumDims)
            KernelDims = kernelDims(this.DistributionNames, NumDims);
            if any(KernelDims)
                % Fill default kernel params:
                if isempty(this.Kernel)
                    this.Kernel = this.DefaultKernel;
                end
                if isempty(this.Support)
                    this.Support = this.DefaultSupport;
                end
                if isempty(this.Width)
                    this.Width = this.DefaultWidth;
                end
                % Check:
                this = checkKernel(this, KernelDims, NumDims);
                this = checkSupport(this, KernelDims, NumDims);
                this = checkWidth(this, KernelDims, NumClasses, NumDims);
            elseif ~isempty(this.Kernel)
                error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:KernelArgButNoKernelDist'));
            elseif ~isempty(this.Support)
                error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:SupportArgButNoKernelDist'));
            elseif ~isempty(this.Width)
                error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:WidthArgButNoKernelDist'));
            end
        end
        
        function this = checkKernel (this, KernelDims, NumDims)
            if ~isempty(this.Kernel)
                if ischar(this.Kernel)
                    [this.Kernel, s] = checkAndCompleteString(this.Kernel, this.LegalKernelNames, 3);
                    if isempty(this.Kernel)
                        error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:KernelUnknown', s));
                    end
                elseif iscell(this.Kernel)
                    % check length
                    if length(this.Kernel) ~= NumDims
                        error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:KernelBadLength'));
                    end
                    % check types
                    if ~all(cellfun(@ischar, this.Kernel(KernelDims)))
                        error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:KernelBadType'));
                    end
                    % check strings
                    for d=1:NumDims
                        if KernelDims(d)
                            [this.Kernel{d}, s] = checkAndCompleteString(this.Kernel{d}, this.LegalKernelNames, 3);
                            if isempty(this.Kernel{d})
                                error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:KernelUnknown', s));
                            end
                        end
                    end
                else
                    error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:KernelBadType'));
                end
            end
        end
        
        function this = checkSupport (this, KernelDims, NumDims)
            if ~isempty(this.Support)
                if ischar(this.Support)
                    [this.Support, s] = checkAndCompleteString(this.Support, this.LegalSupportNames, 3);
                    if isempty(this.Support)
                        error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:SupportUnknown', s));
                    end
                elseif isnumeric(this.Support) && length(this.Support)==2
                    if this.Support(1) >= this.Support(2)
                        error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:SupportBadOrder'));
                    end
                elseif iscell(this.Support)
                    % check length
                    if length(this.Support) ~= NumDims
                        error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:SupportBadLength'));
                    end
                    % check contents of kernel dimensions
                    for d = 1:NumDims
                        if KernelDims(d)
                            if ischar(this.Support{d})   % it's a string
                                [this.Support{d}, s] = checkAndCompleteString(this.Support{d}, this.LegalSupportNames, 3);
                                if isempty(this.Support{d})
                                    error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:SupportUnknown', s));
                                end
                            elseif isnumeric(this.Support{d}) && length(this.Support{d})==2 % it's [L U]
                                if this.Support{d}(1) >= this.Support{d}(2)
                                    error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:SupportBadOrder'));
                                end
                            else % it's neither a string nor [L U] vector
                                error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:SupportEntryBadType'));
                            end
                        end
                    end
                else
                    error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:SupportBadType'));
                end
            end
        end

        function this = checkWidth (this, KernelDims, NumClasses, NumDims)
            if ~isempty(this.Width)
                w = this.Width;
                % check data type
                if ~isnumeric(w)
                    error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:WidthBadType'));
                end
                % Check dimensions and expand into a full matrix to
                % simplify further checking
                if isscalar(w)
                    w = repmat(w, NumClasses, NumDims);
                elseif all(size(w) == [NumClasses 1])
                    w = repmat(w, 1, NumDims);
                elseif all(size(w) == [1 NumDims])
                    w = repmat(w, NumClasses, 1);
                elseif all(size(w) == [NumClasses NumDims])    % It's already a full matrix.
                else
                    error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:WidthBadSize'));
                end
                % check contents: require positive width for all non-NaN
                % kernel dims. We know w is a numeric matrix of the right
                % size at this point.
                kW = w(:,KernelDims);           % widths for kernel dims.
                kW = kW(~isnan(kW));            % non-nan widths for kernel dims.
                if ~all(kW(:) > 0)
                    error(message('stats:classreg:learning:modelparams:NaiveBayesParams:NaiveBayesParams:WidthNotPositive'));
                end
            end
        end
	end
end

function N = numClasses (classSummary)
    N = numel(classSummary.ClassNames);
end

function N = numDims (dataSummary)
    % dataSummary.PredictorNames is either a numeric scalar or a cell
    % array.
    if iscell(dataSummary.PredictorNames) 
        N = numel(dataSummary.PredictorNames);
    else
        N = dataSummary.PredictorNames;
    end
end

function [StringOut, StringIn] = checkAndCompleteString (StringIn, LegalStrings, PrefixLength)
    % Check that StringIn is a member of the cell array LegalStrings, up to
    % the first PrefixLength chars, case-insensitive. If it is, return the
    % full string completion, using the first match. If not, return [].
    % Also return the original string.
    i = find(strncmpi(StringIn, LegalStrings, PrefixLength), 1, 'first');
    if isempty(i)
        StringOut = [];
    else
        StringOut = LegalStrings{i};
    end
end

function KernelDims = kernelDims (DistributionNames, NumDims)
    % Returns a logical vector indicating which dimensions are kernel.
    % Assumes all distribution strings have been completed.
    if ischar(DistributionNames)
        DistributionNames = repmat({DistributionNames}, 1, NumDims);
    end
    KernelDims = strcmpi('kernel',DistributionNames);
end