function [p, se, sse] = levenbergMarquardtKernel( fun, p0, varargin )

    if nargin > 2
        maxIter = varargin{1};
    else
        maxIter = 30;
    end
    
    if nargin > 3
        tolX = varargin{2};
    else
        tolX = .01;
    end
    
    % initialize
    p = p0;
    dp = Inf;
    sse0 = Inf;
    J = []; dy = [];
    lambda = 100/sqrt( p'*p );
        
    iter = 1;
    while any( abs(dp./p0) > tolX ) && iter < maxIter
        
if iter > 1
    sseHat = (dy - J*dp)'*(dy - J*dp);
end
        
        % get new results
        [dy, J] = fun( p );
        sse = dy'*dy;

if iter > 1
    RR = (sse0 - sse)./(sse0 - sseHat);
    disp( ['RR: ' num2str( RR )] )
end
% J = J*4;     
        if sse < sse0
            lambda = .1 * lambda;
%             dp = pinv( J'*J + lambda*diag(diag(J'*J)) )*J'*dy;
            dp = pinv( J'*J + lambda*eye(size(J,2)) )*J'*dy;
        else % reduce last step 
            lambda = 10 * lambda;
            J = J0; dy = dy0; p = p0; sse = sse0;
%             dp = pinv( J'*J + lambda*diag(diag(J'*J)) )*J'*dy;
            dp = pinv( J'*J + lambda*eye(size(J,2)) )*J'*dy;
        end   
% J = J/4;        
        % store last iterations results
        p0 = p; sse0 = sse; J0 = J; dy0 = dy;
        
        % update parameters
        p = p + dp;

        % useful output
        disp( ['SSE: ' num2str(sse) ' Parameters: ' num2str( p0' )] ); 
        
        iter = iter + 1;
    end
    
    % return last parameter
    p = p0; 
    
%     se = sqrt(diag(pinv(J'*J))) * sqrt( dy'*dy/(length(dy)-length(p)) );
    v = (length(dy)-length(p));
    covp = pinv(J'*J) * sse / v ;
    delta = tinv( .975, v ) * sqrt( diag( covp ) );
    ci = [p-delta p+delta];
    
end