function G = pwgranger( Y, Pmax )
    for i = 1:size(Y,2)
        for j = 1:size(Y,2)
            if i~=j
                G(i,j) = nirs.math.grangers(Y(:,i), Y(:,j), Pmax);
            else
                G(i,j) = 0;
            end
        end
    end
    
    G(G<0) = 0;
end

