function [p, sse] = gradientDescentKernel( fun, p0, varargin )

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
    
%     % initialize
%     p = p0;
%     dp = Inf;
%     iter = 1;
% 
%     while any( abs(dp./p) > tolX ) && iter < maxIter 
%         [dy, J] = fun( p );
%         
%         % sum squared errors
%         sse = dy'*dy;
%         
%         % gradient 
%         gradF = -2*J'*dy;
%         
%         % solve to find step: 
%         % F(x+dx) = F(x) + gradF * dx = 0
%         dp = - sse/(gradF'*gradF) * gradF;
% 
%         % useful output
%         disp( ['SSE: ' num2str(sse) ' Parameters: ' num2str( p' )] ); 
%         
%         % update parameter vector
%         p = p + dp;
%         
%         iter = iter + 1;
%     end

    % initialize
    p = p0;
    [dy, J] = fun( p );
    sse = dy'*dy;
    
    iter = 1;
    while iter == 1 || (any( abs(dp./p) > tolX ) && iter < maxIter)
        gradF = -2*J'*dy;
        dp = - sse/(gradF'*gradF) * gradF;
        dy = Inf;
        
        % hard max of 33% change in p
        while any( abs(dp) > abs(p/3) )
            dp = dp/2;
        end
        
        % bisection line search
        dp = dp*2;
        while dy'*dy > sse || all( abs(dp./p) < tolX )
            dp = dp/2;
            [dy, J] = fun( p + dp );
        end
        
        sse = dy'*dy;
        p = p + dp;   
        
        % useful output
        disp( ['SSE: ' num2str(sse) ' Parameters: ' num2str( p' )] ); 
                
        iter = iter + 1;
    end
end

    


