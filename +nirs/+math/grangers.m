function C = grangers( x, y, Pmax )

    n = length(y);
    
    X = nirs.math.lagmatrix(x, 1:Pmax);
    Y = nirs.math.lagmatrix(y, 1:Pmax);
    
    [b1, res1] = nirs.math.stepwise_regression([ones(n, 1) Y], y);
    [b2, res2] = nirs.math.stepwise_regression( ...
        [ones(n, 1) Y(:, 1:length(b1)-1) X], y );

    C = max(log( mad(res1,0)^2 / mad(res2,0)^2 ), 0);
end

