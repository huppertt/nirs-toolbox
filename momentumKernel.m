function [p, sse, covp, ci] = momentumKernel( fun, p0, varargin )

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
    dp = Inf*p;
    sse = Inf;
    sseBest = Inf;
    J = []; dy = [];
    
    v = 0*p;
    
    a = 0.01;
    
    iter = 1;
    while any( abs(dp./p0) > tolX ) && iter < maxIter
       
        p = p + 0.9*v;
        
        % get new results
        [dy, J] = fun( p );
        G = -2*J'*dy;
        
        sse = dy'*dy;   
        
        % store iterations results
        if sse < sseBest;
            pBest = p;
            sseBest = sse;
            JBest = J;
        end
        
         % update
        v = 0.9*v - a*G;
        p = p - a*G;
       
%         p0 = p; sse0 = sse; J0 = J; dy0 = dy; dp0 = dp;
        
        % useful output
        disp( ['SSE: ' num2str(sse) ' Parameters: ' num2str( p' )] ); 
        
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