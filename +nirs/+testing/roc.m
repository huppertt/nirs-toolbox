function [tp,fp,th] = roc( truth, pval )

    [th,I] = sort( pval,'ascend' );
    tp = cumsum( truth(I) ) / sum( truth );
    fp = cumsum( ~truth(I) ) / sum(~truth);

end

