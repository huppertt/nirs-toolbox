function [G, F, df1, df2, p] = mvgc(Y, Pmax)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[yf,f] = nirs.math.innovations(Y,Pmax);

[m, n] = size(Y);
X = []; lst = [];
% model selection
for i = 1:n
    
    % column indices for channel i
    idx = length(lst)+1:length(lst) + length(f{i})-1;
    
    % keep track of them
    lst(idx,1) = i;
    
    % design matrix
    X(:,idx) = nirs.math.lagmatrix(Y(:,i), 1:length(f{i})-1);
end



%Remove an DC shifts in the data
[i,j]=find(X==0);
X=X(max(i)+1:end,:);
m=m-max(i);
Y=Y(max(i)+1:end,:);

    % add a constant term
    X   = [ones(m,1) X];
    lst = [0; lst];
    iX=inv(X'*X)*X';
%     for i = 1:n
%         % unrestricted model (all terms)
%         a   = iX * Y(:,i);
%         ru  = Y(:,i) - X*a;
%         disp(i)
%         for j = 1:n
%             if i == j
%                 % diagonal terms
%                 G(i,j) = 0;
%             else
%                 % channel j terms
%                 l = lst ~= j;
%                 
%                 % restricted model (no j terms)
%                 a  = X(:,l) \ Y(:,i);
%                 rr = Y(:,i) - X(:,l)*a;
%                 
%                 G(i,j)      = log(mad(rr)/mad(ru));
%                 df1(i,j)    = size(X,2) - sum(l);
%                 df2(i,j) 	= size(X,1) - size(X,2);
%                 F(i,j)      = ((mad(rr)^2-mad(ru)^2)/df1(i,j))/(mad(ru)^2/df2(i,j));
%                 F(i,j)      = max(F(i,j),0);
%                 p(i,j)      = fcdf(1/F(i,j), df2(i,j), df1(i,j));
%             end
%         end
%     end

    for j = 1:n
       
       
        % channel j terms
        l = lst ~= j;
                
        iXl=inv(X(:,l)'*X(:,l))*X(:,l)';
        
        for i = 1:n
            if i == j
                % diagonal terms
                G(i,j) = 0;
            else
                % unrestricted model (all terms)
                a   = iX * Y(:,i);
                ru  = Y(:,i) - X*a;

                % restricted model (no j terms)
                a  = iXl * Y(:,i);
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

