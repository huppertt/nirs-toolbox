function [p, sse] = gaussNewtonKernel( fun, p0, varargin )

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
    
    % initialize
    p = p0;
    dp = Inf;
    iter = 1;

    while any( abs(dp./p0) > tolX ) && iter < maxIter
        [dy, J] = fun( p );
        
        % sum squared errors
        sse = dy'*dy;
        
        % gauss newton step
        dp = pinv(J)*dy;
        
        % update parameter vector
        p = p + dp;

        % useful output
        disp( ['SSE: ' num2str(sse) ' Parameters: ' num2str( p' )] ); 
        iter = iter + 1;
    end
    
end