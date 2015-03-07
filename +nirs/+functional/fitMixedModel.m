function S = fitMixedModel( X, Z, y, C )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    nu = size(Z,2);
    nb = size(X,2);
    ny = size(y,1);
    
    if nargin < 4
        C = eye(ny);
    end

    %% initialize
    q = 1; s = 1;
    
    % iterate
    q0 = 1e16; s0 = 1e16;
    while norm([q s] - [q0 s0]) > 1e-6
        q0 = q; s0 = s;
        
        V = q*(Z*Z') + s*C;
        
        L = inv( chol(V) );
        iV = L*L';
        
        b = pinv(X'*iV*X)*X'*iV*y;
        u = q*Z'*iV*(y-X*b);
        
        s = (y-X*b)'*(y-X*b) / (ny-nb);
        q = u'*u / nu;
    end
    
%     tic
%     q = 1; s = 1;
% 
%     nu = size(Z,2);
%     nb = size(X,2);
% 
%     q0 = 1e16; s0 = 1e16;
%     while norm([q s] - [q0 s0]) > 1e-6 && all([s q] > 1e-12)
% 
%         q0 = q; s0 = s;
% 
%         lhs = [     [X'*X/s X'*Z/s];
%                     [Z'*X/s eye(nu)/q+Z'*Z/s]  ];
% 
%         rhs = [1/s*X'*y; 1/s*Z'*y];
% 
%         m = lhs \ rhs;
% 
%         b = m(1:nb);
%         u = m(nb+1:end);
% 
%         s = (y-X*b)'*(y-X*b) / (length(y)-rank(X));
%         q = u'*u / length(u);
% 
%     end
%     toc
    
    S.b = b;
    S.u = u;
    
    S.covb  = pinv(X'*iV*X);
    S.t     = b./sqrt(diag(S.covb));
    S.df    = ny-nb;
    

end

