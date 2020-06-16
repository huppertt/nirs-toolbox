function [stim_vectors,stim_names] = aux2stim(aux,threshold)
% This function fixes some of the stim mark issues on the Cw6 system
if nargin<2 || isempty(threshold)
    threshold = 100;
end
stim_vectors=zeros(size(aux,1),0);
stim_names={};

if isa(aux,'nirs.core.GenericData')
    aux = aux.data;
end

deltaAux = abs(diff(aux));
z=bsxfun(@rdivide,bsxfun(@minus,deltaAux,median(deltaAux)),mad(deltaAux(:)));
if ischar(threshold) && any(strcmpi({'best','max'},threshold))
    threshold = max(z(:));
end
hasmarks=any(z>=threshold);
if any(hasmarks)
    cnt=1; bin=zeros(size(aux,1),0);
    for i=1:size(aux,2)
        if hasmarks(i)
            slocal=aux(:,i);
            slocal=slocal-min(slocal);
            slocal=slocal/max(slocal);           
            bin(:,cnt)=2^(i-1)*(slocal>.5);
            cnt=cnt+1;
            stim_names = [stim_names {sprintf('aux_channel%i',i)}];
        end
    end
    bin=sum(bin,2);
    ub=unique(bin); ub(find(ub==0))=[];
    for i=1:length(ub)
        stim_vectors(find(bin==ub(i)),i)=1;
    end
    stim_vectors(1,:)=0;
    
end

end