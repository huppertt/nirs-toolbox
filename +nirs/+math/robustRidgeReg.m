function S = robustRidgeReg( X, y )

    S = robust_rfit(X,y);
    
end

function w = wfun( r )
    s = mad(r,0) / 0.6745;
    r = r / s / 4.685;
    
    w = (1 - r.^2) .* (abs(r) < 1);
end

function S = robust_rfit( X, y )

    b = X \ y;
    r = y - X*b;
    for i = 1:50
        w = wfun( r );
        
        Xw = bsxfun( @times, w, X );
        yw = w.*y;
        
        S = rfit( Xw, yw );
    end
    
    S.covb  = pinv(X'*X + S.a*eye(length(b))) * (mad(S.r,0)/0.6745)^2;
    
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
%     out.covb = pinv(X'*X) * (mad(out.r,0)/0.6745)^2;
    

end