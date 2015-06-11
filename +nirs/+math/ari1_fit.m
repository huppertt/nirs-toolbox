function [coef, res] = ari1_fit( y, Pmax )
    
    [coef, res] = ar_fit(diff([0; y]), Pmax);
    f = conv([1 -1]', [1; -coef(2:end)]);
    coef = [0; -f(2:end)];
    
end