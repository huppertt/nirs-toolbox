function [mua,mus] = Brain( so2,hbt,lambda )

    c = nirs.getSpectra(lambda);
    
    mua = hbt * so2 * c(:,1) + ...
        hbt * (1-so2) * c(:,2) + ...
        0.7 * c(:,3);
    
    mus = 2.42 * (lambda/500).^-1.611;
    
end

