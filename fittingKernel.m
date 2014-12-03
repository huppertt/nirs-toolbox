function [p, sse, covp, ci] = fittingKernel( fun, p0, varargin )

    if nargin > 2
        maxIter = varargin{1};
    else
        maxIter = 50;
    end
    
    if nargin > 3
        tolX = varargin{2};
    else
        tolX = .001;
    end

    % initialize
    p = p0;
    dp = Inf*p;
    sse = Inf;
    sse0 = sse;
    lambda = .001;
    sseBest = Inf;
    J = []; dy = [];
    
    n = length(p)-2;
    
    iter = 1; count = 1; fail = 0;
    while any( abs(dp./p0) > tolX ) && iter < maxIter
       % get new results
        [dy, J] = fun( p );
        
%         ff = 0;
%         if iter > 1
%             dy = ff*dy0 + (1-ff)*dy;
%             J = ff*J0 + (1-ff)*J;
%         end
%        
%         sse = dy'*dy;
        
        
%         if (max(abs(dp./p0)) < 0.01 || n < ceil(iter/5)) ...
%                 && n < length(dp)
%             n = n + 1;
%         end

%         if iter == 5
%             n = n+1;
%         end
%         
%         if iter == 10
%             n = n+1;
%         end
%         
%         
%         [U,S,V] = svd(J,'econ');
%         U = U(:,1:n);
%         S = diag(S);
%         S = S(1:n);
%         V = V(:,1:n);
%         
%         dp = Inf*p;
%         lambda = 1e-3;
%         while any( abs(dp) > abs(p/5) )
%             lambda = 10 * lambda;
%             
%             dp = V*diag( 1./(S + lambda) )*(U' * dy) / 2;
%         end

        
%         if sse < sse0
%             lambda = lambda/5;
%         else
%             lambda = 2*lambda;
%         end        
        
%         iQ = diag(1./[1 10 1 10 1 1].^2);
%         
%         dp = Inf*p;
%         lambda = lambda/2;
%         lambda = 1e-9;
%         while any( abs(dp) > abs(p/4) )
%             lambda = 5 * lambda;
% %             dp = pinv( J'*J + lambda*iQ )*J'*dy;
% %             dp = pinv( J'*J + lambda*diag(diag(J'*J)) )*J'*dy;
%             dp = pinv( J'*J + lambda*diag(diag(J'*J).*[10 1 10 1 10 10]') )*J'*dy;
% %             dp = pinv( J'*J + lambda*eye(size(J,2)) )* J' * dy;
%         end

% if iter > 0
%     iP = eye(length(dy)); %diag(1./(dy.*dy) ); %* max(10-iter,1) );
%     iQ = diag(1./[2 3 2 3 1 1]).^2;
%     m = [0.7 0.4 0.3 0.2 1.2 1.2]';
%     
%     dp = Inf*p; lambda = lambda/2;
% 	while any( abs(dp) > abs(p/4) )
%         lambda = lambda * 2;
%         dp = pinv(J'*iP*J + lambda*iQ) * (J'*iP*dy + lambda*iQ*(p-m));
%     end
%     
% end
        sse = dy'*dy;
        if sse > sse0 && count < 5
            J = J0; dy = dy0; p = p0; sse = sse0; %dp = dp0;
            
            [dy, J] = fun( p );
            
            dy = dy0 * (count - 1)/count + dy / count;
            J = J0 * (count - 1)/count + J / count;
            count = count + 1;
        elseif sse < sse0
            count = 1;
        else
            fail = 1;
        end
        
        % this parameter alters the diagonals of the regularization term
        % it helps convergence by altering search directions through
        % parameter space; it should not affect the actual solution since
        % we have an over-determined system (# meas > # param)
        q = ones(size(p));
        q(1:2:end) = 2;         % dampen search direction in top layer
        q(2*end/3+1:end) = 10;   % dampen search direction in scattering
        
        % L-curve for choosing lambda
        t = linspace(-10,10,500);
        for i = 1:length(t)
            tmp(:,i) = pinv( J'*J + 10^t(i) * diag( diag(J'*J).*q ) ) * J' * dy;
        end
        [~,idx] = max( abs(diff(log10(sqrt(diag(tmp'*tmp))),2)) );
        dp = tmp(:,idx+2);
        
        % hard max of changing parameters by 50%
        while any( abs(dp./p) > 0.5 )
            idx = idx + 5;
            dp = tmp(:,idx+2);
        end  
        
        % store iterations results
        if sse < sseBest;
            pBest = p;
            sseBest = sse;
            JBest = J;
        end
%         % reset state
%         J = J0; dy = dy0; p = p0; sse = sse0; dp = dp0;
        
        p0 = p; sse0 = sse; J0 = J; dy0 = dy; %dp0 = dp;
        
        % update parameters
        if fail == 1
            dp = 0;
        end
        
        p = p + dp;

        % useful output
        disp( ['SSE: ' num2str(sse) ' Parameters: ' num2str( p0' )] ); 
        disp( ['Next Parameters: ' num2str( p' )] ); 
        
        iter = iter + 1;
    end
    
    % return last parameter
    p = pBest;
    
    %%
    v = (length(dy)-length(p));
    covp = pinv(JBest'*JBest) * sseBest / v ;
    deltaP = tinv( .975, v ) * sqrt( diag( covp ) );
    ci = [p-deltaP p+deltaP]; 
end