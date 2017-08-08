
%   Copyright 2011 The MathWorks, Inc.


classdef DiagonalDiscriminantCalculator < classreg.learning.impl.DiscriminantCalculator
    
    methods
        function this = DiagonalDiscriminantCalculator(mu,invD)
            this = this@classreg.learning.impl.DiscriminantCalculator(mu,invD);
        end
        
        function m = linear(this,X)
            m = bsxfun(@times,X,this.InvD);
        end
        
        function v = quadratic(this,X1,X2)
            v = sum(linear(this,X1).*linear(this,X2),2);
        end
        
        function mah = mahal(this,K,X)
            mah = zeros(size(X,1),numel(K));
            i = 0;
            for k=K
                A = linear(this,bsxfun(@minus,X,this.Mu(k,:)));
                i = i+1;
                mah(:,i) = sum(A.*A,2);
            end
        end
        
        function v = linearCoeffs(this,i,j)
            mu1 = this.Mu(i,:);
            mu2 = this.Mu(j,:);
            v = bsxfun(@times,linear(this,mu1-mu2),this.InvD)';
        end
        
        function sig = sigma(this,d,s,v)
            sig = d.^2;
        end
        
        function invsig = invSigma(this)
            invsig = this.InvD.^2;
        end
        
        function logsig = logDetSigma(this,d,s,v)
            logsig = 2*sum(log(d(d>0)));
        end
    end
    
end
