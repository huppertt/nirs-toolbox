function G = mvgc(Y, Pmax)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    [m, n] = size(Y);
    
    lst = []; 
    X = [];
    for i = 1:n
        a = nirs.math.ar_fit(Y(:,i), Pmax);
        idx = length(lst)+1:length(lst) + length(a)-1;
        lst(idx,1) = i;
        X(:,idx) = nirs.math.lagmatrix(Y(:,i), 1:length(a)-1);
    end
    
    X = [ones(m,1) X];
    lst = [0; lst];
    
    for i = 1:n
        
        %[~, ~, ru] = nirs.math.greedyRegression(X, Y(:,i));
        a = X \ Y(:,i);
        ru = Y(:,i) - X*a;
        for j = 1:n
            if i == j
                G(i,j) = 0;
            else
                l = lst ~= j;
                %[~, ~, rr] = nirs.math.greedyRegression(X(:,l), Y(:,i));
                a = X(:,l) \ Y(:,i);
                rr = Y(:,i) - X(:,l)*a;
                
                G(i,j) = log(mad(rr)/mad(ru));
            end
        end
    end

%     for i = 1:size(Y,2)
%         for j = 1:size(Y,2)
%             if i~=j
%                 G(i,j) = nirs.math.grangers(Y(:,i), Y(:,j), Pmax);
%             else
%                 G(i,j) = 0;
%             end
%             
%             disp([i j])
%         end
%     end
%     
%     G(G<0) = 0;


end

