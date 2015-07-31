function [G, F, df1, df2, p] = mvgc(Y, Pmax)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    [m, n] = size(Y);
    
    X = []; lst = []; 
    
    % model selection
    for i = 1:n
        % fit AR model
        a = nirs.math.ar_fit(Y(:,i), Pmax);
        
        % column indices for channel i
        idx = length(lst)+1:length(lst) + length(a)-1;
        
        % keep track of them
        lst(idx,1) = i;
        
        % design matrix
        X(:,idx) = nirs.math.lagmatrix(Y(:,i), 1:length(a)-1);
    end
    
    % add a constant term
    X   = [ones(m,1) X];
    lst = [0; lst];
    
    for i = 1:n
        % unrestricted model (all terms)
        a   = X \ Y(:,i);
        ru  = Y(:,i) - X*a;
        
        for j = 1:n
            if i == j
                % diagonal terms
                G(i,j) = 0;
            else
                % channel j terms
                l = lst ~= j;
                
                % restricted model (no j terms)
                a  = X(:,l) \ Y(:,i);
                rr = Y(:,i) - X(:,l)*a;
                
                G(i,j)      = log(mad(rr)/mad(ru));
                df1(i,j)    = size(X,2) - sum(l);
                df2(i,j) 	= size(X,1) - size(X,2);
                F(i,j)      = ((mad(rr)^2-mad(ru)^2)/df1(i,j))/(mad(ru)^2/df2(i,j));
                F(i,j)      = max(F(i,j),0);
                p(i,j)      = fcdf(1/F(i,j), df2(i,j), df1(i,j));
            end
        end
    end

end

