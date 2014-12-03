function [p, sse] = BFGSKernel( fun, p0, varargin )

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
    
    % initial iteration
    [dy, J] = fun( p );
    
    % initial hessian
%     B = pinv( J'*J );
    B = eye( size(J,2) );
    
    % sse (i.e., objective function)
    sse = dy'*dy;
    gradSSE = 2*J'*dy;

    % step size
    a = 1;
    
    % parameter update
    dp = a * B*gradSSE;
    p = p + dp;
    
    disp( ['SSE: ' num2str(sse) ' Parameters: ' num2str( p' )] ); 
        
    iter = 1;
    while any( abs(dp./p0) > tolX ) && iter < maxIter
        % store last iterations results
        p0 = p; sse0 = sse; J0 = J; dy0 = dy; B0 = B; gradSSE0 = gradSSE;
        
        % get new results
        [dy, J] = fun( p );

        % check results and update lambda
        sse = dy'*dy;
        if sse < sse0
            
            % update hessian
            gradSSE = 2*J'*dy;
            y = gradSSE - gradSSE0;
            B = B0 + (dp'*y + y'*B*y)*(dp*dp') / (dp'*y)^2 ...
                - (B*y*dp' + dp*y'*B) / (dp'*y);
            
            % update parameters
            a = 1;
            dp = a * B*gradSSE;
            p = p + dp;
            
        else % reduce last step 
            a = a/2;
            p = p0 + a * B*gradSSE;
        
        end   

        % useful output
        disp( ['SSE: ' num2str(sse) ' Parameters: ' num2str( p' )] ); 
        
        iter = iter + 1;
    end
    
end