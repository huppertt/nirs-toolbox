function [coef, res, yhat] = ar_fit( y,Pmax )
    % y(t) = coef(1) + coef(2)*y(t-1) + coef(3)*y(t-2) ... 
    %
    % y is a column vector with the time-series we are fitting
    %
    % Pmax is the max model order (not counting the constant)
    
    n = length(y);
    
    Xf = nirs.math.lagmatrix(y, 1:Pmax);
    Xb = nirs.math.lagmatrix(flipud(y), 1:Pmax);
    
    X = [ones(2*n,1) [Xf; Xb]];
    
    [coef, res] = nirs.math.stepwise_regression(X, [y; flipud(y)]);
    
    yhat = y - res(1:n);
end