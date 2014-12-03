function [p, sse, covp, ci] = trustRegionKernel( fun, p0, varargin )

    if nargin > 2
        maxIter = varargin{1};
    else
        maxIter = 30;
    end
    
    if nargin > 3
        tolX = varargin{2};
    else
        tolX = .001;
    end

    nGradIter = 0;
    
    % initialize
    p = p0;
    dp = Inf*p;
    sse = Inf;
    J = []; dy = [];
    
    % trust region
    %no parameter can change more than this
    delta = p/2; 
    
    iter = 1;
    while any( abs(dp./p0) > tolX ) && iter < maxIter
       % get new results
        [dy, J] = fun( p );
       
        sse = dy'*dy;
        
        if iter > 1
            sseHat = (dy0 - J0*dp0)'*(dy0 - J0*dp0);
            RR = (sse0 - sse)/(sse0 - sseHat);
        else
            RR = .25;
        end
        
        if RR > .5
            % increase trust region
            delta = min( 1.5 * delta, p );
        elseif RR < 0
            % decrease trust region
            delta = 0.75 * delta;
            
            % reset state
            J = J0; dy = dy0; p = p0; sse = sse0; dp = dp0;
        end

        % first iterations are gradient descent
        % this greatly improves convergence
        if iter <= nGradIter
            gradF = -2*J'*dy;
            
            while any(abs(dp) > delta)
                dp = - sse/(gradF'*gradF) * gradF;
                gradF = gradF*2;
            end
                
        else
            % reset lambda
            lambda = .0001;

            dp = Inf*p;
            while any( abs(dp) > delta )
                lambda = 10 * lambda;
                dp = pinv( J'*J + lambda*diag(diag(J'*J)) )*J'*dy;
%             	dp = pinv( J'*J + lambda*eye(size(J,2)) )* J' * dy;
            end
        end
        
        % store iterations results
        p0 = p; sse0 = sse; J0 = J; dy0 = dy; dp0 = dp;
        
        % update parameters
        p = p + dp;

        % useful output
        disp( ['SSE: ' num2str(sse) ' Parameters: ' num2str( p0' )] ); 
        
        iter = iter + 1;
    end
    
    % return last parameter
    p = p0;
    
    %% FIX THIS TO RETURN LAST STATS
    v = (length(dy)-length(p));
    covp = pinv(J0'*J0) * sse0 / v ;
    deltaP = tinv( .975, v ) * sqrt( diag( covp ) );
    ci = [p-deltaP p+deltaP]; 
end