classdef SVMParams < classreg.learning.modelparams.ModelParams
%SVMParams Parameters for binary Support Vector Machine.
%
%   SVMParams properties:
%       Alpha                    - Vector of initial estimates for the
%                                  alpha coefficients.
%       BoxConstraint            - Box constraint.
%       CacheSize                - Cache size, MB.
%       CachingMethod            - Caching algorithm such as 'Queue'.
%       DeltaGradientTolerance   - Tolerance on the gradient difference
%                                  between upper and lower pool violators
%                                  for SMO.
%       GapTolerance             - Tolerance on the feasibility gap.
%       KKTTolerance             - Tolerance on the largest KKT violation.
%       IterationLimit           - Maximal number of iterations.
%       KernelFunction           - One of: 'linear', 'gaussian' (or 'rbf'),
%                                  and 'polynomial'.
%       KernelScale              - Factor by which predictors are divided
%                                  before kernel product for a pair of
%                                  observations is computed.
%       KernelOffset             - Value added to the kernel product for a
%                                  pair of observations.
%       KernelPolynomialOrder    - Order of the polynomial kernel.
%       Nu                       - Nu parameter for one-class learning.
%       NumPrint                 - Print-out period for diagnostic
%                                  messages.
%       OutlierFraction          - Expected fraction of outliers in the
%                                  data.
%       ShrinkagePeriod          - Period for shrinking the active set.
%       Solver                   - One of: 'SMO', 'ISDA', and 'L1QP'.
%       StandardizeData          - Logical flag for standardizing training
%                                  data to zero mean and unit variance.
%       SaveSupportVectors       - True if Alpha and SV's are saved, false if only Beta are saved.
%       VerbosityLevel           - One of: 0, 1 or 2. Increasing the level
%                                  produces diagnostic messages.

%   Copyright 2013-2014 The MathWorks, Inc.
    
    properties
        Alpha = [];
        BoxConstraint = [];
        CacheSize = [];
        CachingMethod = '';
        DeltaGradientTolerance = [];
        GapTolerance = [];
        KKTTolerance = [];
        IterationLimit = [];
        KernelFunction = '';
        KernelScale = [];
        KernelOffset = [];
        KernelPolynomialOrder = [];
        NumPrint = [];
        Nu = [];
        OutlierFraction = [];
        ShrinkagePeriod = [];
        Solver = '';
        StandardizeData = [];
        SaveSupportVectors = [];
        VerbosityLevel = [];
    end
    
    methods(Access=protected)
        function this = SVMParams(alphas,C,nu,...
                cacheSize,cacheAlg,...
                deltagradtol,gaptol,kkttol,...
                outfrac,maxiter,...
                kernelfun,scale,offset,polyorder,...
                dostandardize,solver,shrinkAfter,saveSV,...
                verbose,nprint)
            this = this@classreg.learning.modelparams.ModelParams('SVM','classification');

            this.Alpha                    = alphas;
            this.BoxConstraint            = C;
            this.CacheSize                = cacheSize;
            this.CachingMethod            = cacheAlg;
            this.DeltaGradientTolerance   = deltagradtol;
            this.GapTolerance             = gaptol;
            this.KKTTolerance             = kkttol;
            this.OutlierFraction          = outfrac;
            this.IterationLimit           = maxiter;
            this.KernelFunction           = kernelfun;
            this.KernelScale              = scale;
            this.KernelOffset             = offset;
            this.KernelPolynomialOrder    = polyorder;
            this.Nu                       = nu;
            this.NumPrint                 = nprint;
            this.ShrinkagePeriod          = shrinkAfter;
            this.Solver                   = solver;
            this.StandardizeData          = dostandardize;
            this.SaveSupportVectors       = saveSV;
            this.VerbosityLevel           = verbose;            
        end
    end
    
    methods(Static,Hidden)
        function [holder,extraArgs] = make(type,varargin) %#ok<INUSL>
            % Decode input args
            args = {'alpha' 'boxconstraint' 'nu' ...
                'cachesize' 'cachingmethod' ...
                'kkttolerance' 'gaptolerance' 'deltagradienttolerance' ...
                'outlierfraction' ...
                'iterationlimit' ...
                'kernelfunction' 'kernelscale' 'kerneloffset' 'polynomialorder' ...
                'solver' ...
                'standardize' ...
                'shrinkageperiod' ...
                'savesupportvectors' ...
                'verbose' 'numprint'};
            defs = {     []              []   [] ...
                         []              '' ...
                            []             []                       [] ...
                               [] ...
                              [] ...
                              ''            []             []                [] ...
                      '' ...
                           [] ...
                               [] ...
                                  [] ...
                       []       []};
            [alphas,C,nu,cachesize,cachingmethod,...
                kkttol,gaptol,deltagradtol,outfrac,maxiter,...
                kernelfun,scale,offset,polyorder,...
                solver,dostandardize,shrinkAfter,saveSV,...
                verbose,nprint,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            if ~isempty(alphas) && ...
                    (~isfloat(alphas) || ~isvector(alphas) || ...
                    any(alphas<0) || any(isnan(alphas)) || any(isinf(alphas)))
                error(message('stats:classreg:learning:modelparams:SVMParams:make:BadAlpha'));
            end
            alphas = alphas(:);
            if ~isempty(alphas)
                internal.stats.checkSupportedNumeric('Alpha',alphas,true);
            end
            
            if ~isempty(C) && ...
                    (~isscalar(C) || ~isfloat(C) || C<=0 || isnan(C))
                error(message('stats:classreg:learning:modelparams:SVMParams:make:BadBoxConstraint'));
            end
            
            if ~isempty(nu) && ...
                    (~isscalar(nu) || ~isfloat(nu) || nu<=0 || isnan(nu))
                error(message('stats:classreg:learning:modelparams:SVMParams:make:BadNu'));
            end
            
            if ~isempty(cachingmethod) && ~ischar(cachingmethod)
                error(message('stats:classreg:learning:modelparams:SVMParams:make:BadCachingMethod'));
            end
            
            if ~isempty(cachesize)
                if (~ischar(cachesize) || ~strncmpi(cachesize,'maximal',length(cachesize))) ...
                        && (~isscalar(cachesize) || ~isfloat(cachesize) ...
                        || cachesize<=0  || isnan(cachesize) )
                    error(message('stats:classreg:learning:modelparams:SVMParams:make:BadCacheSize'));
                end
                if ischar(cachesize)
                    cachesize = 'maximal';
                end
            end
            
            if ~isempty(kkttol) && ...
                    (~isscalar(kkttol) || ~isfloat(kkttol) ...
                    || kkttol<0 || isnan(kkttol))
                error(message('stats:classreg:learning:modelparams:SVMParams:make:BadKKTTolerance'));
            end
            
            if ~isempty(gaptol) && ...
                    (~isscalar(gaptol) || ~isfloat(gaptol) ...
                    || gaptol<0 || isnan(gaptol))
                error(message('stats:classreg:learning:modelparams:SVMParams:make:BadGapTolerance'));
            end
            
            if ~isempty(deltagradtol) && ...
                    (~isscalar(deltagradtol) || ~isfloat(deltagradtol) ...
                    || deltagradtol<0 || isnan(deltagradtol))
                error(message('stats:classreg:learning:modelparams:SVMParams:make:BadDeltaGradientTolerance'));
            end
            
            if ~isempty(outfrac) && ...
                    (~isscalar(outfrac) || ~isfloat(outfrac) ...
                    || outfrac<0 || outfrac>=1 || isnan(outfrac))
                error(message('stats:classreg:learning:modelparams:SVMParams:make:BadOutlierFraction'));
            end
            
            if ~isempty(maxiter) && ...
                    (~isscalar(maxiter) || ~isfloat(maxiter) ...
                    || maxiter<=0 || isnan(maxiter))
                error(message('stats:classreg:learning:modelparams:SVMParams:make:BadIterationLimit'));
            end
            
            if ~isempty(kernelfun)
                if ~ischar(kernelfun)
                    error(message('stats:classreg:learning:modelparams:SVMParams:make:KernelFunctionNotString'));
                end
                allowedVals = {'linear' 'gaussian' 'rbf' 'polynomial'};
                tf = strncmpi(kernelfun,allowedVals,length(kernelfun));
                Nfound = sum(tf);
                if     Nfound>1
                    error(message('stats:classreg:learning:modelparams:SVMParams:make:KernelFunctionNotRecognized'));
                elseif Nfound==1
                    if find(tf,1)==3
                        kernelfun = 'gaussian';
                    else
                        kernelfun = allowedVals{tf};
                    end
                else % Nfound==0. Must be user-defined function.
                    if exist(kernelfun,'file')==0
                        error(message('stats:classreg:learning:modelparams:SVMParams:make:KernelFunctionFileNotFound',...
                            kernelfun));
                    end
                end
            end
            
            if ~isempty(scale)
                if (~ischar(scale) || ~strncmpi(scale,'auto',length(scale))) ...
                        && (~isscalar(scale) || ~isfloat(scale) ...
                        || scale<=0  || isnan(scale) )
                    error(message('stats:classreg:learning:modelparams:SVMParams:make:BadKernelScale'));
                end
                if ischar(scale)
                    scale = 'auto';
                end
            end
            
            if ~isempty(offset) && ...
                    (~isscalar(offset) || ~isfloat(offset) ...
                    || offset<0 || isnan(offset))
                error(message('stats:classreg:learning:modelparams:SVMParams:make:BadKernelOffset'));
            end
            
            if ~isempty(polyorder) && ...
                    (~isscalar(polyorder) || ~isfloat(polyorder) ...
                    || polyorder<=0 || isnan(polyorder))
                error(message('stats:classreg:learning:modelparams:SVMParams:make:BadPolynomialOrder'));
            end
            
            if ~isempty(solver)
                if ~ischar(solver)
                    error(message('stats:classreg:learning:modelparams:SVMParams:make:SolverNotString'));
                end
                allowedVals = {'SMO' 'ISDA' 'L1QP'};
                tf = strncmpi(solver,allowedVals,length(solver));
                if sum(tf)~=1
                    error(message('stats:classreg:learning:modelparams:SVMParams:make:SolverNotRecognized'));
                end
                solver = allowedVals{tf};
            end
            
            if ~isempty(dostandardize)
                dostandardize = internal.stats.parseOnOff(dostandardize,'Standardize');
            end
            
            if ~isempty(shrinkAfter) && ...
                    (~isscalar(shrinkAfter) || ~isfloat(shrinkAfter) ...
                    || shrinkAfter<0 || isnan(shrinkAfter))
                error(message('stats:classreg:learning:modelparams:SVMParams:make:BadShrinkagePeriod'));
            end
            
            if ~isempty(saveSV)
                saveSV = internal.stats.parseOnOff(saveSV,'SaveSupportVectors');
            end
            
            if ~isempty(verbose) && ...
                    (~isscalar(verbose) || ~isfloat(verbose) ...
                    || verbose<0 || isnan(verbose))
                error(message('stats:classreg:learning:modelparams:SVMParams:make:BadVerbose'));
            end
            
            if ~isempty(nprint) && ...
                    (~isscalar(nprint) || ~isfloat(nprint) ...
                    || nprint<0 || isnan(nprint))
                error(message('stats:classreg:learning:modelparams:SVMParams:make:BadNumPrint'));
            end
            
            holder = classreg.learning.modelparams.SVMParams(alphas,C,nu,...
                cachesize,cachingmethod,...
                deltagradtol,gaptol,kkttol,...
                outfrac,maxiter,...
                kernelfun,scale,offset,polyorder,...
                dostandardize,solver,shrinkAfter,saveSV,...
                verbose,nprint);
        end
        
        function this = loadobj(obj)
            found = fieldnames(obj);
            
            if ismember('SaveSupportVectors',found) ...
                    && ~isempty(obj.SaveSupportVectors)
                saveSV = obj.SaveSupportVectors;
            else
                saveSV = true;
            end
            
            this = classreg.learning.modelparams.SVMParams(...
                obj.Alpha,obj.BoxConstraint,obj.Nu,...
                obj.CacheSize,obj.CachingMethod,...
                obj.DeltaGradientTolerance,obj.GapTolerance,obj.KKTTolerance,...
                obj.OutlierFraction,obj.IterationLimit,...
                obj.KernelFunction,obj.KernelScale,...
                obj.KernelOffset,obj.KernelPolynomialOrder,...
                obj.StandardizeData,obj.Solver,obj.ShrinkagePeriod,saveSV,...
                obj.VerbosityLevel,obj.NumPrint);
        end
    end
    
    methods(Hidden)
        function this = fillDefaultParams(this,X,Y,W,dataSummary,classSummary) %#ok<INUSL>
            N = size(X,1);
            
            if ~isempty(this.Alpha) && numel(this.Alpha)~=N
                error(message('stats:classreg:learning:modelparams:SVMParams:fillDefaultParams:BadAlpha',N));
            end
            
            if isempty(this.BoxConstraint)
                this.BoxConstraint = 1;
            else
                if numel(classSummary.NonzeroProbClasses)==1
                    if this.BoxConstraint~=1
                        error(message('stats:classreg:learning:modelparams:SVMParams:fillDefaultParams:BadBoxConstraint'));
                    end
                end
            end
            
            % Accept Nu for two-class learning as well. If CROSSVAL is run
            % on a binary SVM model and one class is missing in a fold, SVM
            % fits a one-class model using this value of Nu.
            if isempty(this.Nu)
                this.Nu = 0.5;
            else
                if this.Nu>1
                    error(message('stats:classreg:learning:modelparams:SVMParams:fillDefaultParams:BadNu'));
                end
            end
            
            % If any alphas are passed, make sure they make sense. For
            % one-class learning, require that alphas sum to nu*N, where nu
            % is passed through BoxConstraint and the actual box constraint
            % is always 1. For two-class learning, require that alphas lie
            % below the box constraint. For two-class learning the
            % constraint sum(alphas.*Y)=0 is not enforced.
            if ~isempty(this.Alpha) 
                alphas = this.Alpha;
                if numel(classSummary.NonzeroProbClasses)==1
                    sumAlpha = sum(alphas);
                    
                    if sumAlpha==0
                        error(message('stats:classreg:learning:modelparams:SVMParams:fillDefaultParams:BadAlphaForOneClassLearning',...
                                sprintf('%e',N*this.Nu)));
                    end
                    
                    if abs(sumAlpha - N*this.Nu) > 100*eps(sumAlpha)
                        alphas = alphas*N*this.Nu/sumAlpha;
                        if any(alphas>1)
                            error(message('stats:classreg:learning:modelparams:SVMParams:fillDefaultParams:BadAlphaForOneClassLearning',...
                                sprintf('%e',N*this.Nu)));
                        end
                    end
                else
                    if any(alphas>this.BoxConstraint)
                        maxAlpha = max(alphas);
                        alphas = alphas*this.BoxConstraint/maxAlpha;
                    end
                end
                this.Alpha = alphas;
            end
            
            % Make sure at least one column fits
            if     isa(X,'double')
                oneColSize = ceil(8*N/1024/1024); % MB
            elseif isa(X,'single')
                oneColSize = ceil(4*N/1024/1024); % MB
            end
            
            if     isempty(this.CacheSize)                
                if     isa(X,'double')
                    this.CacheSize = max(1000,oneColSize); % in MB
                elseif isa(X,'single')
                    this.CacheSize = max(1000,oneColSize); % in MB
                end
            elseif strcmpi(this.CacheSize,'maximal') || this.CacheSize==Inf
                if     isa(X,'double')
                    this.CacheSize = ceil(8*N*N/1024/1024);
                elseif isa(X,'single')
                    this.CacheSize = ceil(4*N*N/1024/1024);
                end
            else
                if this.CacheSize < oneColSize
                    error(message('stats:classreg:learning:modelparams:SVMParams:fillDefaultParams:BadCacheSize',...
                        oneColSize));
                end
            end
            
            if isempty(this.CachingMethod)
                this.CachingMethod = 'Queue';
            end
                        
            if     isempty(this.OutlierFraction)
                this.OutlierFraction = 0;
            elseif this.OutlierFraction>0
                if numel(classSummary.NonzeroProbClasses)>1
                    if     isempty(this.Solver)
                        this.Solver = 'ISDA';
                    elseif strcmpi(this.Solver,'L1QP')
                        error(message('stats:classreg:learning:modelparams:SVMParams:fillDefaultParams:QPwithRobustLearning'));
                    end
                end
            end
            
            if isempty(this.Solver)
                this.Solver = 'SMO';
            end
            
            if numel(classSummary.NonzeroProbClasses)==1 ...
                    && strcmpi(this.Solver,'ISDA')
                error(message('stats:classreg:learning:modelparams:SVMParams:fillDefaultParams:ISDAwithOneClassLearning'));
            end
            
            if     strcmpi(this.Solver,'L1QP')
                if isempty(this.KernelOffset)
                    this.KernelOffset = 0;
                end
                if ~isempty(this.KKTTolerance)
                    error(message('stats:classreg:learning:modelparams:SVMParams:fillDefaultParams:BadToleranceParameterForL1QP',...
                        'KKTTolerance'));
                end
                if ~isempty(this.GapTolerance)
                    error(message('stats:classreg:learning:modelparams:SVMParams:fillDefaultParams:BadToleranceParameterForL1QP',...
                        'GapTolerance'));
                end
                if ~isempty(this.DeltaGradientTolerance)
                    error(message('stats:classreg:learning:modelparams:SVMParams:fillDefaultParams:BadToleranceParameterForL1QP',...
                        'DeltaGradientTolerance'));
                end                
            elseif strcmpi(this.Solver,'SMO') 
                if isempty(this.KernelOffset)
                    this.KernelOffset = 0;
                end
                if isempty(this.KKTTolerance)
                    this.KKTTolerance = 0; % not used
                end
                if isempty(this.GapTolerance)
                    this.GapTolerance = 0; % not used
                end
                if isempty(this.DeltaGradientTolerance)
                    this.DeltaGradientTolerance = 0.001;
                end
            elseif strcmpi(this.Solver,'ISDA')
                if isempty(this.KernelOffset)
                    this.KernelOffset = 0.1;
                end
                if isempty(this.KKTTolerance)
                    this.KKTTolerance = 0.001;
                end
                if isempty(this.GapTolerance)
                    this.GapTolerance = 0; % not used
                end
                if isempty(this.DeltaGradientTolerance)
                    this.DeltaGradientTolerance = 0; % not used
                end
            elseif strcmpi(this.Solver,'All2D')
                if isempty(this.KernelOffset)
                    this.KernelOffset = 0.01;
                end
                if isempty(this.KKTTolerance)
                    this.KKTTolerance = 0.001;
                end
                if isempty(this.GapTolerance)
                    this.GapTolerance = 0.001;
                end
                if isempty(this.DeltaGradientTolerance)
                    this.DeltaGradientTolerance = 0.001;
                end
            end
            
            if isempty(this.KernelFunction)
                if numel(classSummary.NonzeroProbClasses)==1
                    this.KernelFunction = 'gaussian';
                    this.SaveSupportVectors = true;
                else
                    this.KernelFunction = 'linear';
                end
            end
            
            if isempty(this.SaveSupportVectors)
                this.SaveSupportVectors = true;
            else
                if ~this.SaveSupportVectors && ~strcmp(this.KernelFunction,'linear')
                    error(message('stats:classreg:learning:modelparams:SVMParams:fillDefaultParams:MustKeepSVforNonlinearKernel'));
                end
            end
            
            if isempty(this.KernelScale)
                this.KernelScale = 1;
            else
                % If user passes his own kernel function, scaling must be
                % handled by that function.
                knownKernels = {'linear' 'gaussian' 'rbf' 'polynomial'};
                isknown = ismember(this.KernelFunction,knownKernels);
                if ~isknown && this.KernelScale~=1
                    error(message('stats:classreg:learning:modelparams:SVMParams:fillDefaultParams:KernelScaleForCustomKernel',...
                        this.KernelFunction));
                end
            end
            
            if isempty(this.KernelPolynomialOrder)
                if strcmpi(this.KernelFunction,'polynomial')
                    this.KernelPolynomialOrder = 3;
                end
            else
                if ~strcmpi(this.KernelFunction,'polynomial')
                    error(message('stats:classreg:learning:modelparams:SVMParams:fillDefaultParams:PolyOrderWithoutPolyKernel'));
                end
            end
            
%             if isempty(this.KernelParameters.SigmoidParameters)
%                 if strcmpi(this.KernelParameters.Function,'sigmoid')
%                     this.KernelParameters.SigmoidParameters = [-1 1];
%                 end                
%             end
            
            if isempty(this.IterationLimit)
                this.IterationLimit = 1e6;
            end
            
            if isempty(this.StandardizeData)
                this.StandardizeData = false;
            end
            
            if isempty(this.ShrinkagePeriod)
                this.ShrinkagePeriod = 0;
            end
            
            if isempty(this.VerbosityLevel)
                this.VerbosityLevel = 0;
            end
            
            if isempty(this.NumPrint)
                if this.VerbosityLevel>0
                    this.NumPrint = 1000;
                else
                    this.NumPrint = 0;
                end
            end
        end
    end    

end
