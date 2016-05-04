function [tp,fp,th] = roc( truth, pval )

    lst=find(truth);
    lstN=find(~truth);
    if(length(lst)~=length(lstN))
        warning('roc has uneven binary profile');
        n=min(length(lst),length(lstN));
        lst=lst(randperm(length(lst),n));
        lstN=lstN(randperm(length(lstN),n));
        pval=pval([lst; lstN]);
        truth=truth([lst; lstN]);
    end
    

    [th,I] = sort( pval,'ascend' );
    tp = cumsum( truth(I) ) / sum( truth );
    fp = cumsum( ~truth(I) ) / sum(~truth);

end

