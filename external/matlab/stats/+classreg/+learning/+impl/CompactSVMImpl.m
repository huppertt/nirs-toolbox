classdef CompactSVMImpl
    
%   Copyright 2013-2014 The MathWorks, Inc.

    properties(SetAccess=protected,GetAccess=public)
        Alpha                       = [];
        Beta                        = [];
        Bias                        = [];
        KernelParameters            = [];
        Mu                          = [];
        NumPredictors               = [];
        Sigma                       = [];
        SupportVectors              = []; 
        SupportVectorLabels         = []; 
    end
   
    methods(Access=protected)
        function this = CompactSVMImpl()
        end
    end
    
    methods
        function f = score(this,X)
            if ~isfloat(X) || ~ismatrix(X)
                error(message('stats:classreg:learning:impl:CompactSVMImpl:score:BadX'));
            end
            internal.stats.checkSupportedNumeric('X',X);
            
            if size(X,2)~=this.NumPredictors
                error(message('stats:classreg:learning:impl:CompactSVMImpl:score:BadXSize',this.NumPredictors));
            end
            
            alphas = this.Alpha;
            betas  = this.Beta;
            bias   = this.Bias;
            
            if isempty(alphas) && isempty(betas)
                f = NaN(size(X,1),1);
                return;
            end
            
            mu = this.Mu;
            if ~isempty(mu) && ~all(mu==0)
                X = bsxfun(@minus,X,mu);
            end
            
            sigma = this.Sigma;
            if ~isempty(sigma) && ~all(sigma==1)
                nonzero = sigma > 0;
                if any(nonzero)
                    X(:,nonzero) = bsxfun(@rdivide,X(:,nonzero),sigma(nonzero));
                end
            end
            
            if isa(X,'double') && isa(this.Bias,'single')
                X = single(X);
            end
            
            if isempty(alphas)
                f = (X/this.KernelParameters.Scale)*betas + this.Bias;
            else
                f = classreg.learning.svmutils.predict(...
                    alphas.*this.SupportVectorLabels,bias,this.SupportVectors,...
                    this.KernelParameters.Function,this.KernelParameters.PolyOrder,...
                    this.KernelParameters.Sigmoid,...
                    this.KernelParameters.Scale,X);
            end
        end
        
        function this = discardSupportVectors(this)
            if ~strcmp(this.KernelParameters.Function,'linear')
                error(message('stats:classreg:learning:impl:CompactSVMImpl:discardSupportVectors:CannotDiscardSVforNonlinearKernel'));
            end
            
            this.Alpha               = [];
            this.SupportVectors      = [];
            this.SupportVectorLabels = [];
        end
    end
    
end
