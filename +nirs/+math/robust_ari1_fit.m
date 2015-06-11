function [coef, res, wres, ymoco] = robust_ari1_fit( y, Pmax )
    yd = diff([0; y]);

    [~, res] = nirs.math.ar_fit(yd, Pmax);
    
    w = wfun(res);
    [coef, ~] = ar_fit(w.*yd, Pmax);
    
    res     = filter([1; -coef(2:end)], 1, yd-coef(1));
    wres    = filter([1; -coef(2:end)], 1, w.*yd-coef(1));
    
    ymoco   = filter(1, [1; -coef(2:end)], w.*res);
    ymoco   = cumsum(ymoco);
    
    f = conv([1 -1]', [1; -coef(2:end)]);
    coef = [0; -f(2:end)];
end

function w = wfun(r)
    s = mad(r, 0) / 0.6745;
    r = r/s/4.685;
    
    w = (1 - r.^2) .* (r < 1 & r > -1);
end