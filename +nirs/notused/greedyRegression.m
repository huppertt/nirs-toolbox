function [b, idx, r] = greedyRegression(X, y)
    [m, n] = size(X);
    
    [Q, R] = qr(X,0);
    
    iR = pinv(R);
    
    lst = ones(n,1) > 0;
    idx = zeros(n,1);
    
    for i = 1:n
        b = iR(lst, lst) * Q(:,lst)' * y;
        r = y - X(:,lst)*b;
        
        covb = iR(lst,lst)*iR(lst,lst)' * var(r);
        t = b ./ sqrt(diag(covb));
        
        % calculate information criterion
        LL = -m/2*log( 2*pi*var(r) ) - m/2; % log-likelihood
        k = n-i+1;
        bic(i) = -LL + k/2*log(m);
        %aic(i) = -LL + k + k*(k+1)/(m-k-1);
        
        % remove least significant variable
        [~, j] = min(t.^2);
        ii = find(lst);
        j = ii(j);
        
        lst(j) = false;
        idx(i) = j;
    end
    
    [~, k] = min(bic);
    
    idx = sort(idx(k:end));
    
    b = iR(idx, idx) * Q(:,idx)' * y;
    r = y - X(:,idx)*b;
end