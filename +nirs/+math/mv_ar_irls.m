function S = mv_ar_irls( X, Y, Pmax )

    [m, n, p] = size(X);
    
    stack = @(X) reshape( permute(X,[1 3 2]), [m*p n]);
    vec   = @(Y) Y(:);
    
    % initial fit
    S = nirs.math.robust_mvfit( stack(X), Y );
    
    % residual
    r = vec(Y) - stack(X)*S.b;
    r = reshape(r, size(Y));
    
    % whiten data and model
    Yf = zeros(size(Y));
    Xf = zeros(size(X));
    for i = 1:size(Y,2)
        a           = ar_fit( r(:,i), Pmax );
        Yf(:,i)     = myFilter( [1; a(2:end)], Y(:,i) );
        Xf(:,:,i)   = myFilter( [1; a(2:end)], X(:,:,i) );
    end
        
    % regress
    S = nirs.math.robust_mvfit( stack(Xf), Yf );

end

function out = myFilter( f, y )
    % here we are just making the first value zero before filtering to
    % avoid weird effects introduced by zero padding
    y1 = y(1,:);
    
    y = bsxfun(@minus,y,y1);
    
    out = filter(f, 1, y);
    out = bsxfun(@plus,out,sum(f)*y1); % add the corrected offset back

end