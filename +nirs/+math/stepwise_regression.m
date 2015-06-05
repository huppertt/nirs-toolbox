function [b, r] = stepwise_regression( X, y)

 	% qr factorization will speed up stepwise regression significantly
    [Q,R] = qr(X,0); % zero is very important for performance
    invR = pinv(R);  % save this for permance reasons
    
    BIC = zeros(size(X,2),1);
    for i = 1:length(BIC)
        % get residual for each fit for orders 0 to Pmax
        b = invR(1:i,1:i) * Q(:,1:i)'*y;
        r = y - X(:,1:i)*b;
        
        % calculate information criterion
        n = length(r);
        LL = -n/2*log( 2*pi*var(r) ) - n/2; % log-likelihood
        
        % uncomment the line you want to use for information criterion
        % I find that BIC produces more consistent results, but some
        % studies suggest AIC is superior
        BIC(i) = -LL + i/2*log(n);
    end
    
    x = linspace(-1, 1, length(BIC))';
    x = [x.^0 x.^1 x.^2];
    BIC = x*(x\BIC);
    
    [~, N] = min( BIC ); % optimal model order
    
    % finally, our output
    b = invR(1:N,1:N) * Q(:,1:N)'*y;
    r = y - X(:,1:N)*b;
    
end

