function tbdata=get_stimdesign_aux(taux,aux);

tbdata={};
cnt=1;
for idx=1:size(aux,2)
    [onsets,durations]=getstim(taux,aux(:,idx));
    for id=1:length(onsets)
    tbdata{cnt,1}=onsets(id);
    tbdata{cnt,2}=durations(id);
    tbdata{cnt,3}=1;
    tbdata{cnt,4}=['aux_' num2str(idx)];
    
    cnt=cnt+1;
    end
    
    
end



return

function [Onsets,Durations]=getstim(taux,aux);

thres=.25;
minT=4;
fs=1/mean(diff(taux));

aux=aux-mean(aux);
aux=zscore(aux);
if(mean(aux(find(taux<3)))>thres)
    dl=find(diff(aux)<-thres/2);
    aux(1:min(dl))=0;
end

Onsets=find(diff(aux)>thres);
Onsets(find(diff(Onsets)<minT*fs))=[];

Onsets(find(Onsets>(length(aux)-minT*fs)))=[];

Durations=[];
for idx=1:length(Onsets)
    oo=find(diff(aux)<-thres);
    oo(find(oo<Onsets(idx)))=[];
    Durations(idx)=min(oo); 
    
end

if(~isempty(Onsets))
    Onsets=taux(Onsets)';
    Durations=taux(Durations)'-Onsets;
end

return