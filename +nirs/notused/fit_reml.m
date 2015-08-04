function out = fit_reml( X, y, Sc, Pc )
    if nargin < 4
        Pc = {eye(size(X,2))};
    end
    
    if nargin < 3
        Sc = {eye(size(X,1))};
    end
    
    opts = optimset('MaxIter',3);
    
    % fit model
    s = ones( length(Sc), 1 );
    p = ones( length(Pc), 1 );
    
    for i = 1:3
        [S, P]  = getCovMats( Sc, Pc, s, p );
        
        % fit b given q, r
        [b, covb] = remlInv( X, y, S, P );
        
        % fit q given b, r
        c = @(s) nll( X, y, b, Sc, Pc, s, p );
        s = fmincon(c, s, [], [], [], [], 0*s, Inf*s);
        
        % fit r given b, q
        c = @(p) nll( X, y, b, Sc, Pc, s, p );
        p = fmincon(c, p, [], [], [], [], 0*p, Inf*p);
        
        disp(i)
    end
    
    out.b       = b;
    out.covb    = covb;
    out.t       = b ./ sqrt( diag( covb ) );
    out.s       = s;
    out.p       = p;
    
end

function out = nll( X, y, b, Sc, Pc, s, p )
    
    [S, P] = getCovMats( Sc, Pc, s, p );
    
    iS = S \ eye(size(X,1));
    iP = P \ eye(size(X,2));
    
   	out = (y-X*b)' * iS * (y-X*b) ...
        + b' * iP * b ...
        + log( det(S) ) ...
        + log( det(P) );
    
end

function [S, P] = getCovMats( Sc, Pc, s, p )
    S = 0; P = 0;
    for i = 1:length(s)
        S = S + s(i)*Sc{i};
        P = P + p(i)*Pc{i};
    end
end

function [b, C] = remlInv( X, y, S, P )
    
    iS = S \ eye(size(X,1));
    iP = P \ eye(size(X,2));
    
    C = pinv( X'*iS*X + iP );
    b =  C * X' * iS * y;
    
end
