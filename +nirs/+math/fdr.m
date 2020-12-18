function q = fdr( p )
% See http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530395/
    s = size(p);
    
    % sort the p values
    [p, I] = sort(p(:));

    % number of hypothesis tests
    m = length(p);

    if m < 100
        % Benjamini-Hochberg Procedure
        pi0 = 1;
        
    elseif m < 1000
        % Storey's method w/ point estimate
        pi0 = 0.5 / ( 1 - p(round(end/2)) );
        
    else
        % Storey's method w/ interpolation
        % grid to estimate fraction of null tests
        x = p;
        y = ( m - (1:m)' ) / m ./ (1-x);
        
        % only in this range
        lst = x > 0.1 & x < 0.85;
        y = y(lst);
        x = x(lst);
    
        % estimate
        X = [ones(size(x)) x x.^2 x.^3];
        b = X \ y;

        % interpolate
        pi0 = [1 1 1 1] * b;
 
        pi0 = max(pi0, 0);
        pi0 = min(pi0, 1);
    end
    pi0=1;
     % p = (i/m)*Q
    q = p*m./[1:m]'*pi0;
        
    for i=length(q):-1:2; if(q(i)<q(i-1)); q(i-1)=q(i); end; end;
    
    % put them back in the original order
    q(I) = q;
    
    q = reshape(q, s);
    
    q( q > 1 ) = 1;
end

