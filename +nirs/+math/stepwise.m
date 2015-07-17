function [b, r, bic] = stepwise(X, y)
    % qr factorization will speed up stepwise regression significantly
    [Q,R] = qr(X,0); % note that the zero is very important for performance
    invR = pinv(R);
    
    n = length(y);
    bic = zeros(size(X,2),1);
    for i = 1:length(bic)
        % get residual for each fit
        b = invR(1:i,1:i) * Q(:,1:i)' * y;
        r = y - X(:,1:i)*b;
        
        % calculate information criterion
        LL = -n/2*log( 2*pi*var(r) ) - n/2; % log-likelihood
        bic(i) = -LL + i/2*log(n);
    end
    
    % optimal model order
    [~, N] = min( bic ); 
    
    % finally, our output
    b = invR(1:N,1:N) * Q(:,1:N)'*y;
    r = y - X(:,1:N)*b;
    
end

