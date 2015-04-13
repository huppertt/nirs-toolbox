function L = legendre( t, P )

    assert( isvector(t) && isscalar(P) )
    
    % polynomials will be centered
    t = linspace(-1,1,length(t))';
    
    % preallocation
    L = zeros(size(t,1),P+1);
    
    % 0th order
    L(:,1) = 1; 
    
    if P > 0
        L(:,2) = t; % 1st order

        for i = 2:P % 2nd to Pth order
            L(:,i+1) = ( (2*i+1) * L(:,i) .* t - i * L(:,i-1) ) / (i+1);
        end
    end

end

