function S = fitMixedModel( X, Z, y, C )
% This performs maximum likelihood estimation of a mixed effects model
% specified by X and Z.  C is the covariance of y.

    nu = size(Z,2); % num rfx
    nb = size(X,2); % num ffx
    ny = size(y,1); % num data pts
    
    if nargin < 4
        C = eye(ny);
    end

    %% initialize
    q = 1; % variance of random effects
    s = 1; % variance of errors
    
    %% iterate
    q0 = 1e16; s0 = 1e16;
    while norm([q s] - [q0 s0]) > 1e-6
        q0 = q; s0 = s;
        
        % covariance of y
        V = q*(Z*Z') + s*C;
        
        % invert V
        L = inv( chol(V) );
        iV = L*L';
        
        % fixed effects
        b = pinv(X'*iV*X)*X'*iV*y;
        
        % random effects
        u = q*Z'*iV*(y-X*b);
        
        % update s and q
        s = (y-X*b)'*(y-X*b) / (ny-nb);
        q = u'*u / nu;
    end
    
    %% output
    S.b = b;
    S.u = u;
    
    S.covb  = pinv(X'*iV*X);
    S.t     = b./sqrt(diag(S.covb));
    S.df    = ny-nb;
    

end

