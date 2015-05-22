function q = fdr( p )
% See http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530395/

    % sort the p values
    [p, I] = sort(p);

    % number of hypothesis tests
    m = length(p);
    
    if m > 300 % Storey's method for large number of tests
        
        % grid to estimate fraction of null tests
        x = (0.05:0.05:0.95)';

        y = zeros(size(x));
        for i = 1:length(x)
           y(i) = sum( p > x(i) ) / m / (1-x(i)); 
        end

        % interpolate to y( x = 1 )
        % this is the estimated fraction of null samples
        X = [ones(size(x)) x x.^2 x.^3];
        b = X \ y;

        pi0 = [1 1 1 1] * b;
        pi0 = max(pi0, 0);
        pi0 = min(pi0, 1);

        % this is the q-value, defined as the minimum FDR cutoff
        % that would reject this test
        q = p * m ./ (1:m)' * pi0;
        
    else % BH: not enough samples to estimate pi0
        
        q = p * m ./ (1:m)';
        
    end
    
    % put them back in the original order
    q(I) = q;
end

