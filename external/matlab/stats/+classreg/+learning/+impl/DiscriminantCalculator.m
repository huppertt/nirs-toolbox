
%   Copyright 2011 The MathWorks, Inc.


classdef DiscriminantCalculator < classreg.learning.internal.DisallowVectorOps
    
    properties(GetAccess=public,SetAccess=public)
        % Class means, K-by-p
        Mu = [];
        
        % Inverse of the standard deviation per predictor, 1-by-p
        InvD = [];
    end

    methods(Access=protected)
        function this = DiscriminantCalculator(mu,invD)
            this = this@classreg.learning.internal.DisallowVectorOps();
            this.Mu = mu;
            this.InvD = invD(:)';
        end
    end
    
    methods(Abstract)
        % X/D/Corr, N-by-p. D is the diagonal matrix with standard
        % deviations per predictor and Corr is the correlation matrix.
        % (Covariance matrix is given by Sigma = D*Corr*D)
        m = linear(this,X)
        
        % X1*inv(Sigma)*X2', N-by-1
        v = quadratic(this,X1,X2)
        
        % Square of Mahalanobis distance for class indices in K,
        % N-by-numel(K)
        m = mahal(this,K,X)
        
        % Linear coefficients for classes i and j, p-by-1
        v = linearCoeffs(this,i,j)
        
        % Covariance matrix
        sig = sigma(this,d,s,v)
        
        % Inverse of the covariance matrix
        invsig = invSigma(this)
        
        % Log of the determinant of the covariance matrix
        logsig = logDetSigma(this,d,s,v)
    end
    
    methods
        % Const term for separation of classes i and j
        function c = constantTerm(this,i,j)
            mah = mahal(this,[i j],zeros(1,size(this.Mu,2)));
            c = 0.5*(mah(1)-mah(2));
        end
    end
end
