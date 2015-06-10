function L = legendre( t, P )

    assert( isvector(t) && isscalar(P) )

    % recursive definition: not orthogonal
    % gram-schmidt: not numerically stable
    
    %% using cholesky decomposition
    n = length(t);
    
    x = linspace(-1, 1, n)';                        % polys on [-1 1]
    
    L = bsxfun( @power, x, 0:P );                   % polynomials of x
    L = bsxfun( @rdivide, L, sqrt(sum(L.^2,1)) );   % normalize
    
    L = L*inv(chol(L'*L));                          % orthonormalization
    

end