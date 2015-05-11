function S = robustRidgeReg( X, y )

    % initial fit
    S = rfit( X, y );

    % loop
    count = 0; b0 = 1e16;
    while norm(S.b-b0)/norm(b0) > 1e-4 && count < 50
        b0 = S.b;
        
        w = wfun( S.r );
        
        Xw = bsxfun( @times, w, X );
        yw = w.*y;
        
        S = rfit( Xw, yw );
        
        count = count + 1;
    end
    
    S.covb  = pinv(X'*X + 0*S.a*eye(length(S.b))) * (mad(S.r,0)/0.6745)^2;
    S.dfe = size(X,1) - size(X,2);
    S.w = w;
end

function w = wfun( r )
    s = mad(r,0) / 0.6745;
    r = r / s / 4.685;
    
    w = (1 - r.^2) .* (abs(r) < 1);
end

function out = rfit( X, y )

    a = 10.^(-6:6);
    
    [U, S, V] = svd(X, 'econ');
    s = diag(S);
    
    for i = 1:length(a)
        b = V * diag(s ./ (s.^2 + a(i)^2)) * U' * y;
        
        r = y - X*b;
        df = size(X,1) - sum( s ./ (s.^2 + a(i)^2) );
        
        gcv(i) = r'*r / df^2;
    end
    
    [~,i] = min(gcv);
    
    out.b = V * diag(s ./ (s.^2 + a(i)^2)) * U' * y;
    out.r = y - X*out.b;
    out.a = a(i);    

end