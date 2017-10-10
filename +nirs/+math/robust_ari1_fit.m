function [coef, res, wres, ymoco] = robust_ari1_fit( y, Pmax, tune )
    if nargin < 3, tune = 4.685; end
    
    yd = diff([0; y]);

    [~, res] = nirs.math.ar_fit(yd, Pmax);
    
    w = wfun(res, tune);
    [coef, ~] = nirs.math.ar_fit(w.*yd, Pmax);
    
    res     = filter([1; -coef(2:end)], 1, yd-coef(1));
    wres    = w.*res; %filter([1; -coef(2:end)], 1, w.*yd-coef(1));
    
    f = conv([1 -1]', [1; -coef(2:end)]);
    
    yf = filter(f,1,y);
    w = wfun(yf, tune);
    ymoco = filter(1,f,w.*yf);
    
    coef = [0; -f(2:end)];
end

function w = wfun(r, tune)
if(1)
    r=nirs.math.normrootstationarity(r,'mean');
    r=nirs.math.normrootstationarity(r,'std');
    r=r-nanmean(r);
    lstN=find(r<0);
    lstP=find(r>0);
    
    s = mad(r(lstN), 0) / 0.6745;
    r(lstN) = r(lstN)/s/tune;
    
    s = mad(r(lstP), 0) / 0.6745;
    r(lstP) = r(lstP)/s/tune;
else
    %original version
    s = mad(r, 0) / 0.6745;
    r = r/s/tune;
end
w = (1 - r.^2) .* (r < 1 & r > -1);
w(isnan(w))=0;
end