classdef SVMImpl < classreg.learning.impl.CompactSVMImpl
    
%   Copyright 2013-2014 The MathWorks, Inc.

    properties(SetAccess=protected,GetAccess=public)
        Active                      = [];
        C                           = [];
        CacheInfo                   = [];
        ConvergenceInfo             = [];
        FractionToExclude           = [];
        Gradient                    = [];
        IsSupportVector             = [];
        NumIterations               = 0;
        NumMatrixVectorProductTaken = 0;
        PreferredUpdate             = [];
        Shrinkage                   = [];
        WorkingSet                  = [];
    end
   
    methods(Access=protected)
        function this = SVMImpl()
            this = this@classreg.learning.impl.CompactSVMImpl();
        end
    end
    
    methods
        function cmp = compact(this,saveSV)
            cmp = classreg.learning.impl.CompactSVMImpl();
            cmp.Beta                = this.Beta;
            cmp.Bias                = this.Bias;
            cmp.KernelParameters    = this.KernelParameters;
            cmp.Mu                  = this.Mu;
            cmp.NumPredictors       = this.NumPredictors;
            cmp.Sigma               = this.Sigma;            
            if saveSV
                cmp.Alpha               = this.Alpha;
                cmp.SupportVectors      = this.SupportVectors;
                cmp.SupportVectorLabels = this.SupportVectorLabels;
            end
        end
        
        function this = resume(this,X,Y,numIter,doclass,verbose,nprint)
        %RESUME Resume training by SMO or ISDA
        %   OBJ=RESUME(OBJ,X,Y,NUMITER,DOCLASS,VERBOSE,NPRINT) resumes training the
        %   model OBJ on data X and Y for NUMITER iterations. Pass Y as a vector
        %   filled with -1 and +1. Pass DOCLASS as true or false. Pass VERBOSE as a
        %   non-negative integer. Pass NPRINT as a non-negative integer.

            if isa(X,'single')
                Y = single(Y);
            end
            
            if this.PreferredUpdate==0
                % The model must've been obtained by quadprog
                error(message('stats:classreg:learning:impl:SVMImpl:resume:CannotResumeForSolver'));
            end
            
            if this.ConvergenceInfo.Converged
                error(message('stats:classreg:learning:impl:SVMImpl:resume:CannotResumeAfterConvergence'));
            end
            
            mu = this.Mu;
            if ~isempty(mu) && ~all(mu==0)
                X = bsxfun(@minus,X,mu);
            end
            sigma = this.Sigma;
            if ~isempty(sigma) && ~all(sigma==1)
                nonZeroSigma = sigma>0;
                if any(nonZeroSigma)
                    X(:,nonZeroSigma) = ...
                        bsxfun(@rdivide,X(:,nonZeroSigma),sigma(nonZeroSigma));
                end
            end
            
            N = size(X,1);
            
            alphas = zeros(N,1);
            alphas(this.IsSupportVector) = this.Alpha;
            
            if isa(X,'double')
                alphaTol = sqrt(eps);
            else % single
                alphaTol = 10*sqrt(eps);
            end
            
            grad = this.Gradient;
            active = this.Active;
            
            kernelFun     = this.KernelParameters.Function;
            polyOrder     = this.KernelParameters.PolyOrder;
            sigmoidParams = this.KernelParameters.Sigmoid;
            kernelScale   = this.KernelParameters.Scale;
            kernelOffset  = this.KernelParameters.Offset;
            
            wssAlgs = this.WorkingSet.Algorithms;
            preferredUpdate = this.PreferredUpdate;

            c = this.C;
            
            maxIter = this.NumIterations + numIter;
            nFirstIter = this.NumIterations;
            
            kktTol             = this.ConvergenceInfo.KKTTolerance;
            gapTol             = this.ConvergenceInfo.GapTolerance;
            deltaGradTol       = this.ConvergenceInfo.DeltaGradientTolerance;
            outlierHistory     = this.ConvergenceInfo.OutlierHistory;
            
            nIterHistory       = this.ConvergenceInfo.History.NumIterations;
            gapHistory         = this.ConvergenceInfo.History.Gap;
            deltaGradHistory   = this.ConvergenceInfo.History.DeltaGradient;
            worstViolHistory   = this.ConvergenceInfo.History.LargestKKTViolation;
            nsvHistory         = this.ConvergenceInfo.History.NumSupportVectors;
            objHistory         = this.ConvergenceInfo.History.Objective;
            changeSetHistory   = this.ConvergenceInfo.ChangeSetHistory;

            cacheSize = this.CacheInfo.Size;
            cachingAlg = this.CacheInfo.Algorithm;

            nShrinkAfter = this.Shrinkage.Period;
            shrinkAlgs   = this.Shrinkage.Algorithms;
            
            fExclude = this.FractionToExclude;
                        
            [alphas,active,grad,bias,nIter,...
                reasonForConvergence,...
                gap,deltaGradient,kktViolation,Q,...
                nMtimesv,~,wssCounts,...
                outlierHistory,...
                nIterHistory,gapHistory,deltaGradHistory,...
                worstViolHistory,nsvHistory,objHistory,changeSetHistory] = ...
                classreg.learning.svmutils.solve(...
                alphas,grad,active,...
                X,Y,...
                kernelFun,polyOrder,sigmoidParams,kernelScale,kernelOffset,...
                wssAlgs,preferredUpdate,...
                c,maxIter,alphaTol,kktTol,gapTol,deltaGradTol,...
                cacheSize,cachingAlg,nShrinkAfter,shrinkAlgs,fExclude,...
                nFirstIter,outlierHistory,...
                nIterHistory,gapHistory,deltaGradHistory,...
                worstViolHistory,nsvHistory,objHistory,changeSetHistory,...
                verbose,nprint);
            
            if doclass==1
                bias = -1 - quantile(grad,fExclude);
                Q = Q + sum(alphas);
            end
            
            idxSV = alphas>0;

            ws.Algorithms = this.WorkingSet.Algorithms;
            ws.Names      = this.WorkingSet.Names;
            if     isempty(wssCounts)
                ws.Counts = this.WorkingSet.Counts;
            elseif isempty(this.WorkingSet.Counts)
                ws.Counts = wssCounts;
            else
                ws.Counts = this.WorkingSet.Counts + wssCounts;
            end
            
            history.NumIterations       = nIterHistory;
            history.Gap                 = gapHistory;
            history.DeltaGradient       = deltaGradHistory;
            history.LargestKKTViolation = worstViolHistory;
            history.NumSupportVectors   = nsvHistory;
            history.Objective           = objHistory;

            this.Active                               = active;
            this.Alpha                                = alphas(idxSV);
            this.Bias                                 = bias;
            this.ConvergenceInfo.Converged            = ~strcmpi(reasonForConvergence,'NoConvergence');
            this.ConvergenceInfo.ReasonForConvergence = reasonForConvergence;
            this.ConvergenceInfo.Gap                  = gap;
            this.ConvergenceInfo.DeltaGradient        = deltaGradient;
            this.ConvergenceInfo.LargestKKTViolation  = kktViolation;
            this.ConvergenceInfo.Objective            = Q;
            this.ConvergenceInfo.OutlierHistory       = outlierHistory;
            this.ConvergenceInfo.History              = history;
            this.ConvergenceInfo.ChangeSetHistory     = changeSetHistory;
            this.Gradient                             = grad;
            this.IsSupportVector                      = idxSV;
            this.NumIterations                        = nIter;
            this.NumMatrixVectorProductTaken          = this.NumMatrixVectorProductTaken + nMtimesv;
            this.WorkingSet                           = ws;
            this.SupportVectors                       = X(idxSV,:);
            this.SupportVectorLabels                  = Y(idxSV);
        end
    end
    
    methods(Static)
        function this = make(X,Y,W,alphas,...
                kernelFun,polyOrder,sigmoidParams,kernelScale,kernelOffset,...
                doscale,doclass,solver,c,nu,...
                maxIter,kktTol,gapTol,deltaGradTol,...
                cacheSize,cachingAlg,...
                nShrinkAfter,fExclude,verbose,nprint)
            
            % doclass:
            %    0 - regression
            %    1 - one-class classification
            %    2 - two-class classification

            if isa(X,'single')
                Y = single(Y);
                W = single(W);
                c = single(c);
            end
                                    
            if ~strcmp(kernelFun,'polynomial')
                polyOrder = 0; % dummy just to pass in something
            end
            
             if ~strcmp(kernelFun,'builtin-sigmoid') 
                sigmoidParams = [0 0]; % dummy just to pass in something
            end            
            
            [N,D] = size(X);
            
            mu    = [];
            sigma = [];
            if doscale
                mu = classreg.learning.internal.wnanmean(X,W);
                sigma = classreg.learning.internal.wnanvar(X,W,1);
                X = bsxfun(@minus,X,mu);
                nonZeroSigma = sigma>0;
                sigma(~nonZeroSigma) = 0;
                if any(nonZeroSigma)
                    sigma = sqrt(sigma);
                    X(:,nonZeroSigma) = ...
                        bsxfun(@rdivide,X(:,nonZeroSigma),sigma(nonZeroSigma));
                end
            end
            
            if     strcmpi(solver,'SMO')
                wssAlgs = {'MaxViolatingPair' 'MaxDeltaQ'};
                shrinkAlgs = {'UpDownGradient'};
                preferredUpdate = 2;
            elseif strcmpi(solver,'ISDA')
                wssAlgs = {'WorstViolator' 'MaxGainFromHalves'};
                shrinkAlgs = {'KKTViolators'};
                preferredUpdate = 1;
            elseif strcmpi(solver,'All2D')
                wssAlgs = {'MaxDeltaQ' ... 
                    'MaxGainAndPrevFound' 'MaxGainFromHalves' 'MaxGainAndNearby'};
                shrinkAlgs = {'KKTViolators'};
                preferredUpdate = 2;
            else 
                preferredUpdate = 0;
            end

            if strcmpi(kernelScale,'auto')
                kernelScale = classreg.learning.svmutils.optimalKernelScale(...
                    X,Y,doclass);
            end

            kernelParams.Function  = kernelFun;
            kernelParams.PolyOrder = polyOrder;
            kernelParams.Sigmoid   = sigmoidParams;
            kernelParams.Scale     = kernelScale;
            kernelParams.Offset    = kernelOffset;
            
            active = true(N,1);
            grad   = -ones(N,1);
            
%             if     doclass==0
%                 c = repmat(c,N,1);
%             else
            if doclass==2
                c = W*N*c;
            end

            alphasAtBC = false;
            if isempty(alphas)    
                if doclass==1
                    alphas = repmat(nu,N,1); % sum(alphas)=nu*N
                    if nu==1
                        alphasAtBC = true;
                    end
                else
                    alphas = zeros(N,1);
                end
            end
            
            if any(alphas>0)
                if ~strcmpi(solver,'L1QP') || alphasAtBC
                    % Compute the initial gradient
                    idxSV = alphas>0;
                    y_times_alphas = Y(idxSV).*alphas(idxSV);
                    grad = Y.*classreg.learning.svmutils.predict(...
                        y_times_alphas,kernelOffset*sum(y_times_alphas),X(idxSV,:),...
                        kernelFun,polyOrder,sigmoidParams,kernelScale,X) - 1;
                end
            end
            
            if doclass==1
                c = ones(N,1);
                fExcludeOneClass = fExclude;
                fExclude = 0; % do not remove observations with large gradients
            end
            
            % Tolerance on alpha coefficients.
            % Experimentally, these values work best.
            if isa(X,'double')
                alphaTol = sqrt(eps);
            else % single
                alphaTol = 10*sqrt(eps);
            end

            if     alphasAtBC
                % Handle special case: all alphas are at the box constraint
                % for one-class learning
                convergenceInfo.Converged                  = true;
                convergenceInfo.ReasonForConvergence       = ...
                    getString(message('stats:classreg:learning:impl:SVMImpl:make:AllAlphasMustBeAtBoxConstraint'));
                convergenceInfo.GapTolerance               = gapTol;
                convergenceInfo.DeltaGradientTolerance     = deltaGradTol;
                convergenceInfo.KKTTolerance               = kktTol;
                nIter = 0;
                nMtimesv = 0;
                shrinkage = [];
                ws = [];
                Q = alphas'*(grad-1)/2;
                
            elseif strcmpi(solver,'L1QP')
                % Check Optim license
                if isempty(ver('Optim'))
                    error(message('stats:classreg:learning:impl:SVMImpl:make:QPNeedsOptim'));
                end
                
                cachingAlg = '';
                nMtimesv = [];
                shrinkage = [];
                ws = [];
                
                G = classreg.learning.svmutils.computeGramMatrix(...
                    X,kernelFun,polyOrder,sigmoidParams,kernelScale,kernelOffset);
                G = (G+G')/2;
                
                % Increase tolerance for post-QP processing
                oldAlphaTol = alphaTol;
                alphaTol = 100*alphaTol;
                
                if any(c<alphaTol)
                    idxbad = find(c<alphaTol,1);
                    error(message('stats:classreg:learning:impl:SVMImpl:make:QPBadBoxConstraintForObservation',...
                        idxbad,sprintf('%e',c(idxbad))));
                end
                
                opts = optimoptions(@quadprog,...
                    'Algorithm','interior-point-convex',...
                    'TolX',oldAlphaTol,'TolCon',oldAlphaTol);
                opts.MaxIter = maxIter;
                if     verbose>0
                    opts.Display = 'iter';
                else
                    opts.Display = 'none';
                end
                
                if     doclass==1
                    [alphas,Q,exitflag,output] = ...
                        quadprog(double(G),[],[],[],...
                        ones(1,N),double(sum(alphas)),...
                        zeros(N,1),double(c),[],opts);
                elseif doclass==2
                    [alphas,Q,exitflag,output] = ...
                        quadprog(double((Y*Y').*G),-ones(N,1),[],[],...
                        double(Y)',double(sum(alphas.*Y)),...
                        zeros(N,1),double(c),[],opts);
                end
                
                if exitflag <= 0
                    error(message('stats:classreg:learning:impl:SVMImpl:make:NoQPConvergence',...
                        sprintf('\n %s',output.message)));
                end
                
                convergenceInfo.Converged      = exitflag==1;
                convergenceInfo.QuadprogOutput = output;
                
                % Find at least 2 SV's
                while sum(alphas>alphaTol) < 2
                    alphaTol = alphaTol/10;
                end
                idxSV = alphas>alphaTol;
                
                nIter = output.iterations;
                
                y_times_alphas = Y(idxSV).*alphas(idxSV);
                grad = Y.*classreg.learning.svmutils.predict(...
                    y_times_alphas,kernelOffset*sum(y_times_alphas),X(idxSV,:),...
                    kernelFun,polyOrder,sigmoidParams,kernelScale,X) - 1;
                
                % Are there free support vectors?
                isFree = alphas>alphaTol & alphas<c-alphaTol;
                
                % If there are free SV's, average over them. Otherwise
                % average over max gradients for the two groups of
                % violators.
                if any(isFree)
                    bias = mean(-Y(isFree).*grad(isFree));
                else
                    isUp   = (Y==1 & alphas<=c-alphaTol) | (Y==-1 & alphas>=alphaTol);
                    isDown = (Y==1 & alphas>=alphaTol)   | (Y==-1 & alphas<=c-alphaTol);
                    
                    maxUp   = max(-Y(isUp).*grad(isUp));
                    minDown = min(-Y(isDown).*grad(isDown));
                    
                    bias = (maxUp+minDown)/2;
                end
                
            else % SMO or ISDA, regular case
                if any(c<alphaTol)
                    idxbad = find(c<alphaTol,1);
                    error(message('stats:classreg:learning:impl:SVMImpl:make:BadBoxConstraintForObservation',...
                        idxbad,sprintf('%e',c(idxbad))));
                end
                
                nFirstIter = 0;
                [alphas,active,grad,bias,nIter,...
                    reasonForConvergence,...
                    gap,deltaGradient,kktViolation,Q,...
                    nMtimesv,wssNames,wssCounts,...
                    outlierHistory,...
                    nIterHistory,gapHistory,deltaGradHistory,...
                    worstViolHistory,nsvHistory,objHistory,changeSetHistory] = ...
                    classreg.learning.svmutils.solve(...
                    alphas,grad,active,...
                    X,Y,kernelFun,polyOrder,sigmoidParams,kernelScale,kernelOffset,...
                    wssAlgs,preferredUpdate,...
                    c,maxIter,alphaTol,kktTol,gapTol,deltaGradTol,...
                    cacheSize,cachingAlg,nShrinkAfter,shrinkAlgs,fExclude,...
                    nFirstIter,[],...
                    [],[],[],[],[],[],[],... % histories
                    verbose,nprint);
                
                idxSV = alphas>0;
                
                convergenceInfo.Converged                  = ~strcmpi(reasonForConvergence,'NoConvergence');
                convergenceInfo.ReasonForConvergence       = reasonForConvergence;
                convergenceInfo.Gap                        = gap;
                convergenceInfo.GapTolerance               = gapTol;
                convergenceInfo.DeltaGradient              = deltaGradient;
                convergenceInfo.DeltaGradientTolerance     = deltaGradTol;
                convergenceInfo.LargestKKTViolation        = kktViolation;
                convergenceInfo.KKTTolerance               = kktTol;
                convergenceInfo.OutlierHistory             = outlierHistory;
                convergenceInfo.History.NumIterations      = nIterHistory;
                convergenceInfo.History.Gap                = gapHistory;
                convergenceInfo.History.DeltaGradient      = deltaGradHistory;
                convergenceInfo.History.LargestKKTViolation= worstViolHistory;
                convergenceInfo.History.NumSupportVectors  = nsvHistory;
                convergenceInfo.History.Objective          = objHistory;
                convergenceInfo.ChangeSetHistory           = changeSetHistory;
                
                ws.Algorithms = wssAlgs;
                ws.Names      = wssNames;
                ws.Counts     = wssCounts;
                
                shrinkage.Period     = nShrinkAfter;
                shrinkage.Algorithms = shrinkAlgs;
            end
            
            if doclass==1
                % f = grad + bias + 1;
                % Find bias such that the fExcludeOneClass fraction of the
                % one class has negative scores.
                bias = -1 - quantile(grad,fExcludeOneClass);
                
                % Correct the objective by the sum of the alpha
                % coefficients (constant).
                Q = Q + sum(alphas);
                
                % Copy the excluded fraction back to fExclude for storage
                fExclude = fExcludeOneClass;
            end
            
            convergenceInfo.Objective = Q;
            
            beta = [];
            if strcmpi(kernelFun,'linear')
                if any(idxSV)
                    beta = sum( bsxfun( @times, X(idxSV,:), ...
                        Y(idxSV).*alphas(idxSV)), 1 )' / kernelScale;
                    % beta = X(idxSV,:)'*(Y(idxSV).*alphas(idxSV)); % simpler but needs transposition
                else
                    beta = [];
                end
            end
            
            this = classreg.learning.impl.SVMImpl();

            this.Active                       = active;
            this.Alpha                        = alphas(idxSV);
            this.Beta                         = beta;
            this.Bias                         = bias;
            this.C                            = c;
            this.CacheInfo.Size               = cacheSize;
            this.CacheInfo.Algorithm          = cachingAlg;
            this.ConvergenceInfo              = convergenceInfo;
            this.FractionToExclude            = fExclude;
            this.Gradient                     = grad;
            this.IsSupportVector              = idxSV;
            this.KernelParameters             = kernelParams;
            this.Mu                           = mu;
            this.NumIterations                = nIter;
            this.NumMatrixVectorProductTaken  = nMtimesv;
            this.NumPredictors                = D;            
            this.PreferredUpdate              = preferredUpdate;
            this.Shrinkage                    = shrinkage;
            this.Sigma                        = sigma;
            this.WorkingSet                   = ws;
            this.SupportVectors               = X(idxSV,:);
            this.SupportVectorLabels          = Y(idxSV);
        end
    end
    
end
