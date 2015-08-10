function S = mv_ar_irls( X, Y, Pmax, T )

    if nargin < 4, T = ones(size(Y,1),1); end
    
    if ndims(X) == 3
        [m, n, p] = size(X);
    else
        [m, p] = size(Y);
        n = size(X,2);
    end

    unstack = @(X) permute( reshape(X, [m n p]), [1 3 2] );
    stack	= @(X) reshape( permute(X, [1 3 2]), [m*p n] );
    vec     = @(Y) Y(:);
    
    if ndims(X) == 2
        X = unstack(X);
    end
    
    % initial fit
    b = mySolve([stack(X) kron(eye(size(Y,2)),T)], vec(Y));
            
    b0 = b * 1e16; iter = 0;
    while norm(b-b0)/norm(b0) > 1e-3 && iter < 10
        
%         disp( norm(b-b0)/norm(b0) )
        
        b0 = b;
        
        % residual
        r = vec(Y) - [stack(X) kron(eye(size(Y,2)),T)]*b;
        r = reshape(r, size(Y));

        % whiten data and model
        Yf = zeros(size(Y));
        Xf = zeros(size(X));
        Tf = zeros([m*p p*size(T,2)]);
        
        for i = 1:size(Y,2)
            a = nirs.math.ar_fit( r(:,i), Pmax );
            f = [1; -a(2:end)];
       
            Yf(:,i)     = myFilter( f, Y(:,i) );
            Xf(:,i,:)   = myFilter( f, squeeze(X(:,i,:)) );
            
            idx1 = (i-1)*m+1 : i*m;
            idx2 = (i-1)*size(T,2)+1 : i*size(T,2);
            
            Tf(idx1, idx2) = myFilter( f, T );
        end
        
        % filtered residual
        rf = reshape( vec(Yf) - [stack(Xf) Tf]*b, size(Yf) );

        % weight by covariance of resid
        [u, s, ~] = svd(cov(rf),'econ');

        % this is a whitening filter when multiplied on the right
        % i.e. cov(r * W) = I
        Q = u * diag( 1./sqrt(diag(s)) );
%         Q = diag(1./std(rf,1));
%         Q = inv(chol(cov(rf0)));

        % spatial prewhitening
        Yq = Yf*Q;
        Xq = myKronProd( Q, stack(Xf) );
        Tq = myKronProd( Q, Tf );
        
        for j = 1:10
            % tukey weights
            rq = reshape( vec(Yq) - [Xq Tq] * b, size(Yq) );

            w  = wfun( rq );

            % weighted fit
            Xw = bsxfun(@times, Xq, w(:));
            Tw = bsxfun(@times, Tq, w(:));

            Yw = w.*Yq;

            b = mySolve([Xw Tw], vec(Yw));
        end
        
        % update iter count
        iter = iter + 1;
    end
    
    S.b     = b;
    S.covb  = ([Xw Tw]'*[Xw Tw]) \ eye(size(b,1));
    S.dfe   = numel(Y) - numel(b);
end

function b = mySolve( X, y )
    b = (X'*X) \ (X'*y);
end

function X = myKronProd( Q, X )
    % calculates  kron( Q', I ) * X   more efficiently
    vec = @(x) x(:);
    
    for i = 1:size(X,2)
        X(:,i) = vec( reshape(X(:,i), [size(X,1)/size(Q,1) size(Q,1)]) * Q );
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

function w = wfun( x )
    c = 4.685;
    s = mad(x(:), 0) / 0.6745;
    w = (1 - (x/s/c).^2) .* (abs(x/s) < c);
end