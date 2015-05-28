function S = mv_ar_irls( X, Y, Pmax, T )
    if nargin < 4, T = ones(size(Y,1),1); end

    [m, n, p] = size(X);
    
    stack = @(X) reshape( permute(X,[1 3 2]), [m*p n]);
    vec   = @(Y) Y(:);
    
    % initial fit
    thisX = [stack(X) kron(eye(size(Y,2)),T)];
    S = nirs.math.robust_mvfit( thisX, Y );
    
    b0 = S.b * 1e16; iter = 0;
    while norm(S.b-b0) / norm(b0) > 1e-4 && iter < 50
        b0 = S.b;
        
        % residual
        r = vec(Y) - thisX*S.b;
        r = reshape(r, size(Y));

        % whiten data and model
        Yf = zeros(size(Y));
        Xf = zeros(size(X));
        Tf = [];
        for i = 1:size(Y,2)
            a           = ar_fit( r(:,i), Pmax );
            f           = [1; -a(2:end)];
            Yf(:,i)     = myFilter( f, Y(:,i) );
            Xf(:,:,i)   = myFilter( f, X(:,:,i) );
            Tf          = blkdiag(Tf, myFilter( f, T ));
        end

        % regress
        S = nirs.math.robust_mvfit( [stack(Xf) Tf], Yf );
        
        iter = iter + 1;
    end

end

function out = myFilter( f, y )
    % here we are just making the first value zero before filtering to
    % avoid weird effects introduced by zero padding
    y1 = y(1,:);
    
    y = bsxfun(@minus,y,y1);
    
    out = filter(f, 1, y);
    out = bsxfun(@plus,out,sum(f)*y1); % add the corrected offset back

end