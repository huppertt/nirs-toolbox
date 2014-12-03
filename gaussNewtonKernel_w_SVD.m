function [p, sse] = gaussNewtonKernel_w_SVD( fun, p0, varargin )

    if nargin > 2
        maxIter = varargin{1};
    else
        maxIter = 10;
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
        
        [U,S,V] = svd(J,'econ');
        
        S = 1./diag(S);
        
        S = S(1:end-1);
        U = U(:,1:end-1);
        V = V(:,1:end-1);

        % gauss newton step
        dp = (bsxfun(@times,V,S.')*U')*dy;
        
        % update parameter vector
        p = p + dp;

        % useful output
        disp( ['SSE: ' num2str(sse) ' Parameters: ' num2str( p' )] ); 
        iter = iter + 1;
    end
    
end