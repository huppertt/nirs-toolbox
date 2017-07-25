function [b, r, crit] = stepwise(X, y, criterion)
    
    if nargin<3, criterion = 'BIC'; end
    
    % qr factorization will speed up stepwise regression significantly
    [Q,R] = qr(X,0); % note that the zero is very important for performance
    invR = pinv(R);
    
    n = size(y,1);
    crit = zeros(size(X,2),1);
    for i = 1:length(crit)
        % get residual for each fit
        b = invR(1:i,1:i) * Q(:,1:i)' * y;
        r = y - X(:,1:i)*b;
        
        % calculate information criterion
        LL = -n/2*log( 2*pi*mean(r.^2) ) - n/2; % log-likelihood
        
        switch upper(criterion)
            case 'BIC'
                crit(i) = -2*LL + i*log(n);                  % BIC
            case 'AIC'
                crit(i) = -2*LL + 2*i;                       % AIC
            case 'AICC'
                crit(i) = -2*LL + 2*i + 2*i*(i+1)/(n-i-1);   % AICc
            case 'CAIC'
                crit(i) = -2*LL + i*(1+log(n));              % CAIC
            case 'MAX'
                crit(i) = -i;
            otherwise
                error('Unknown model selection criterion: %s',criterion);
        end
    end
    
    % optimal model order
    [~, N] = min( crit ); 
    
    % finally, our output
    b = invR(1:N,1:N) * Q(:,1:N)'*y;
    r = y - X(:,1:N)*b;
    
end

