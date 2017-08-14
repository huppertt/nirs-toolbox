function s = aux2stim(data)
% This function fixes some of the stim mark issues on the Cw6 system
s=data.s;

if(~isfield(data,'aux'))
    return
end

deltaAux = abs(diff(data.aux));
hasmarks=any(bsxfun(@rdivide,bsxfun(@minus,deltaAux,median(deltaAux)),mad(deltaAux))>1000);
if any(hasmarks)
    cnt=1; bin=zeros(length(data.t),0);
    ss=zeros(length(data.t),0);
    for i=1:size(data.aux,2)
        if hasmarks(i)
            slocal=data.aux(:,i);
            slocal=slocal-min(slocal);
            slocal=slocal/max(slocal);           
            bin(:,cnt)=2^(i-1)*(slocal>.5);
            cnt=cnt+1;
        end
    end
    bin=sum(bin,2);
    ub=unique(bin); ub(find(ub==0))=[];
    for i=1:length(ub)
        ss(find(bin==ub(i)),i)=1;
    end
    ss(1,:)=0;
     s=ss;   
else
    s=data.s;
end

return