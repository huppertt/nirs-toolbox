
%   Copyright 2011 The MathWorks, Inc.


classdef DiscriminantImpl < classreg.learning.internal.DisallowVectorOps
    
    properties(GetAccess=public,SetAccess=protected)
        % Type of this discriminant. For LDA can be 'linear', 'diagLinear'
        % or 'pseudoLinear'.
        Type = '';
        
        % K-by-P matrix of class means for K classes and P predictors.
        Mu = [];
        
        % Vector of class weights with K elements.
        ClassWeights = [];
        
        % Global mean vector, 1-by-P
        BetweenMu = [];
        
        % Mu-BetweenMu, K-by-P
        CenteredMu = [];
        
        % Regularization parameters
        Gamma = [];
        Delta = [];
    end
    
    properties(GetAccess=public,SetAccess=public,Abstract=true,Hidden=true)
        SaveMemory;
    end
    
    properties(GetAccess=public,SetAccess=protected,Abstract=true,Hidden=true)
        Sigma;
        InvSigma;
        LogDetSigma;
        MinGamma;
    end
    
    properties(Abstract,Constant)
        % The derived class should set this property to a cell array with
        % all acceptable values. For LDA, for instance, this would be
        % {'linear' 'diagLinear' 'pseudoLinear'}
        AllowedTypes;
    end
    
    methods(Abstract)
        m = linear(this,X)
        v = quadratic(this,X1,X2)
        m = mahal(this,K,X)
        v = linearCoeffs(this,i,j)
        c = constantTerm(this,i,j)

        delran = deltaRange(this,gamma)
        delpred = deltaPredictor(this,gamma)
        
        nCoeffs = nLinearCoeffs(this,delta)
        
        this = setType(this,type)
        this = setGamma(this,gamma)
        this = setDelta(this,delta)
    end
    
    methods
        function this = DiscriminantImpl(gamma,delta,mu,classWeights)
            this = this@classreg.learning.internal.DisallowVectorOps();
            this.Mu = mu;
            this.Gamma = gamma;
            this.Delta = delta;
            [betweenMu,centeredMu,classWeights] = ...
                classreg.learning.impl.DiscriminantImpl.centerMu(mu,classWeights);
            this.ClassWeights = classWeights;
            this.BetweenMu = betweenMu;
            this.CenteredMu = centeredMu;
        end
    end
            
    methods(Static)
        function [betweenMu,centeredMu,Wj] = centerMu(mu,classWeights)
        % Return the global mean vector, class means centered by the global
        % mean and class weights with weights for zero-probability classes
        % set to NaN.
            
            % Classes with zero probabilities have NaN means.
            % Assume that a class has either an entire row filled with
            % NaN's or an entire row of valid numbers.
            Wj = classWeights;
            tfnan = isnan(mu);
            tfall = all(tfnan,2);
            tfany = any(tfnan,2);
            if any(tfall~=tfany)
                error(message('stats:classreg:learning:impl:DiscriminantImpl:centerMu:BadMu'));
            end
            Wj(tfall) = NaN;
            Wj(~tfall) = Wj(~tfall)/sum(Wj(~tfall));
            betweenMu = sum(bsxfun(@times,mu(~tfall,:),Wj(~tfall)),1);
            centeredMu = bsxfun(@minus,mu,betweenMu);
        end
    end
    
end
