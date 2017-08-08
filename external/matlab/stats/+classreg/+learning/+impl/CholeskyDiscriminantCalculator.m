
%   Copyright 2011 The MathWorks, Inc.


classdef CholeskyDiscriminantCalculator < classreg.learning.impl.DiscriminantCalculator
    
    properties(GetAccess=public,SetAccess=public)
        % Inverse of the R, where Corr = R'*R
        InvR = [];
    end
    
    methods
        function this = CholeskyDiscriminantCalculator(mu,invD,invR)
            this = this@classreg.learning.impl.DiscriminantCalculator(mu,invD);
            this.InvR = invR;
        end
        
        function m = linear(this,X)
            m = (bsxfun(@times,X,this.InvD)*this.InvR)*this.InvR';
        end
        
        function v = quadratic(this,X1,X2)
            m1 = bsxfun(@times,X1,this.InvD)*this.InvR;
            m2 = bsxfun(@times,X2,this.InvD)*this.InvR;
            v = sum(m1.*m2,2);
        end
        
        function mah = mahal(this,K,X)
            mah = zeros(size(X,1),numel(K));
            i = 0;
            for k=K
                A = bsxfun(@times,bsxfun(@minus,X,this.Mu(k,:)),this.InvD)*this.InvR;
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
            R = bsxfun(@times,s,bsxfun(@times,v',d));
            sig = R'*R;
        end
        
        function invsig = invSigma(this)
            invR = bsxfun(@times,this.InvD',this.InvR);
            invsig = invR*invR';
        end
        
        function logsig = logDetSigma(this,d,s,v)
            logsig = 2*sum(log(d(d>0))) + 2*sum(log(abs(s(abs(s)>0))));
        end
    end
    
end
