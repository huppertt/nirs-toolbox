function [p, sse, covp, ci] = fittingKernel( fun, p0, varargin )

    if nargin > 2
        maxIter = varargin{1};
    else
        maxIter = 100;
    end
    
    if nargin > 3
        tolX = varargin{2};
    else
        tolX = .0001;
    end

    % initialize
    p = p0;
    dp = Inf*p;
    sse = Inf;
    sse0 = sse;
    lambda = 0.001;
    sseBest = Inf;
    J = []; dy = [];
    
    n = 3;
    
    iter = 1;
    while any( abs(dp./p0) > tolX ) && iter < maxIter
       % get new results
        [dy, J] = fun( p );
        
        ff = 0.9;
        if iter > 1
            dy = ff*dy0 + (1-ff)*dy;
            J = ff*J0 + (1-ff)*J;
        end
       
        sse = dy'*dy;
        
        
        if (max(abs(dp./p0)) < 0.02 || n < ceil(iter/10)) ...
                && n < length(dp)
            n = n + 1;
        end
        
        
        [U,S,V] = svd(J,'econ');
        U = U(:,1:n);
        S = diag(S);
        S = S(1:n);
        V = V(:,1:n);
        
        dp = Inf*p;
        lambda = 1e-3;
        while any( abs(dp) > abs(p/3) )
            lambda = 10 * lambda;
            
            dp = V*diag( 1./(S + lambda) )*(U' * dy) / 10;
        end

        
%         if sse < sse0
%             lambda = lambda/2;
%         else
%             lambda = 2*lambda;
%         end        
%         
%         dp = Inf*p;lambda = lambda/10;
%         while any( abs(dp) > abs(p/4) )
%             lambda = 10 * lambda;
%             dp = pinv( J'*J + lambda*diag(diag(J'*J)) )*J'*dy;
% %             	dp = pinv( J'*J + lambda*eye(size(J,2)) )* J' * dy;
%         end
        
        % store iterations results
        if sse < sseBest;
            pBest = p;
            sseBest = sse;
            JBest = J;
        end
%         % reset state
%         J = J0; dy = dy0; p = p0; sse = sse0; dp = dp0;
        
        p0 = p; sse0 = sse; J0 = J; dy0 = dy; dp0 = dp;
        
        % update parameters
        p = p + dp;

        % useful output
        disp( ['SSE: ' num2str(sse) ' Parameters: ' num2str( p0' )] ); 
        
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