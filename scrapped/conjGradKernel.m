function [p, sse] = conjGradKernel( fun, p0, varargin )

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
    sse = Inf;
    J = []; dy = []; gradF = [];
%     lambda = 1/sqrt( p'*p );
        
    iter = 1;
    while any( abs(dp./p0) > tolX ) && iter < maxIter
        % store last iterations results
        p0 = p; sse0 = sse; J0 = J; dy0 = dy; dp0 = dp; gradF0 = gradF;
        
        % get new results
        [dy, J] = fun( p );
        
        % check results and update lambda
        sse = dy'*dy;
%         if sse < sse0
%             lambda = .1 * lambda;
%             dp = pinv( J'*J + lambda*diag(diag(J'*J)) )*J'*dy;
%         else % reduce last step 
%             lambda = 10 * lambda;
%             J = J0; dy = dy0; p = p0; sse = sse0;
%             dp = pinv( J'*J + lambda*diag(diag(J'*J)) )*J'*dy;
%         end  

%         if sse < sse0
if iter > 1
            gradF = -2*J'*dy;
            beta = (gradF'*gradF)/(gradF0'*gradF0);
            
            newDir = gradF0 + beta*gradF;
            alpha = mean( pinv( newDir/-2*pinv(dy) )'*dy );
            p = p + alpha*newDir;
else
            gradF = -2*J'*dy;
            dp = pinv(J)*dy;
            alpha = mean(dp./gradF);
            p = p + dp;
end

%             dp = pinv( J'*J + lambda*diag(diag(J'*J)) )*J'*dy;
%         else % reduce last step 
%             J = J0; dy = dy0; p = p0; sse = sse0;
% %             dp = pinv( J'*J + lambda*diag(diag(J'*J)) )*J'*dy;
%         end 
        
        % useful output
        disp( ['SSE: ' num2str(sse) ' Parameters: ' num2str( p' )] ); 
        
        iter = iter + 1;
    end
    
end