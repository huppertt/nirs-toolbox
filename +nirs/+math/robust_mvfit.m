function S = robust_mvfit( X, Y )
    [m, n] = size(Y);

    vec = @(x) x(:);
    
    % initial estimate
    b = X \ vec(Y);
    
    b0 = 1e16*b; iter = 0;
    while norm(b - b0)/norm(b0) > 1e-3 && iter < 50
        b0 = b;
        
        % residual & covariance
        r = reshape( vec(Y) - X*b, [m n] );
        R  = cov(r);

        % weight by covariance of resid
        [U,s,~] = svd(R,'econ');
        
        % this is a whitening filter when multiplied on the right
        % i.e. cov(r * W) = I
        W = U * diag( 1./sqrt(diag(s)) );

        Yw = Y*W;                       % whiten on the right
        Xw = kron(W', speye(m)) * X;    % on the left

        % call robustfit to get the stats
        [b, sts] = robustfit(Xw, vec(Yw), [], [], 'off');

        iter = iter + 1;
    end
    
    S.b    = b;
    S.covb = sts.covb;
    S.dfe  = sts.dfe;
end

