
%   Copyright 2011 The MathWorks, Inc.


classdef RegularizedDiscriminantCalculator < classreg.learning.impl.DiscriminantCalculator
    
    properties(GetAccess=public,SetAccess=public)
        % Inverse of the correlation matrix
        InvCorr = [];        

        % Ridge parameter
        Gamma = [];
        
        % Lasso parameter
        Delta = [];
        
        % Global class mean, a 1-by-p vector
        BetweenMu = [];
        
        % Class means with BetweenMu subtracted and divided by the standard
        % deviation per predictor, K-by-p
        CenteredScaledMu = [];
        
        % CenteredScaledMu/Corr, K-by-p
        CenteredScaledMuOverCorr = [];
    end
    
    methods
        function this = RegularizedDiscriminantCalculator(...
                mu,invD,invCorr,gamma,delta,betweenMu,centeredMu)
            this = this@classreg.learning.impl.DiscriminantCalculator(mu,invD);
            this.InvCorr = invCorr;
            this.Gamma = gamma;
            this.Delta = delta;
            this.BetweenMu = betweenMu;            
            this.CenteredScaledMu = bsxfun(@times,centeredMu,this.InvD);
            this.CenteredScaledMuOverCorr = this.CenteredScaledMu*this.InvCorr;
        end
        
        function m = linear(this,X)
            m = bsxfun(@times,X,this.InvD)*this.InvCorr;
        end
        
        function v = quadratic(this,X1,X2)
            v = sum(linear(this,X1).*bsxfun(@times,X2,this.InvD),2);            
        end
        
        function mah = mahal(this,K,X)
            mah = zeros(size(X,1),numel(K));
            standardX = bsxfun(@times,bsxfun(@minus,X,this.BetweenMu),this.InvD);
            qX = sum((standardX*this.InvCorr).*standardX,2);            
            i = 0;
            for k=K
                mu = this.CenteredScaledMuOverCorr(k,:);
                mu(abs(mu)<this.Delta) = 0;
                i = i+1;
                mah(:,i) = qX - bsxfun(@minus,2*standardX,this.CenteredScaledMu(k,:))*mu';
            end
        end
        
        function v = linearCoeffs(this,i,j)
            Rm1 = this.CenteredScaledMuOverCorr(i,:);
            Rm2 = this.CenteredScaledMuOverCorr(j,:);
            Rm1(abs(Rm1)<this.Delta) = 0;
            Rm2(abs(Rm2)<this.Delta) = 0;
            v = bsxfun(@times,Rm1-Rm2,this.InvD)';
        end
        
        function sig = sigma(this,d,s,v)
            gamma = this.Gamma;
            s = sqrt((1-gamma)*s.^2 + gamma);
            R = bsxfun(@times,s,bsxfun(@times,v',d));
            sig = R'*R;
        end
        
        function invsig = invSigma(this)
            invsig = bsxfun(@times,this.InvD',bsxfun(@times,this.InvCorr,this.InvD));
        end
        
        function logsig = logDetSigma(this,d,s,v)
            gamma = this.Gamma;
            logsig = 2*sum(log(d)) + sum(log((1-gamma)*s.^2+gamma));
        end
    end
    
end
