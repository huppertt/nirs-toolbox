function data = loadBrainVision(filename,fs)
% This function loads BrainVision EEG data

if(nargin<2)
    fs=250;
end

[hdr.fs, hdr.label,hdr.meta]=bva_readheader(filename);

n=round(hdr.fs/fs);

d=bva_loadeeg(filename);
aux=d(33:end,:);

d=resample(double(d'),1,n)';
fs=hdr.fs/n;

data=eeg.core.Data;
data.data=d(1:31,:)';
data.time=[0:size(d,2)-1]/fs;
data.description=which(filename);

data.probe=eeg.core.Probe({hdr.label{1:31}});

if(exist([strtok(filename,'.') '.mat']))
    load([strtok(filename,'.') '.mat']);
    stim=nirs.util.convertStimDesignStruct(StimDesign);
else
    stim=findstim(aux,hdr.fs,data.time);
end
data.stimulus=stim;
data.description=filename;
end

function stim = findstim(aux,fs,t)

aux2=aux;
stim=Dictionary;
for i=1:size(aux,1)
    s=diff(aux(i,:));
    s=s-s(1);
    s=s./sqrt(var(s));
    lst=find(s>20);
    lst(find(diff(lst)<50))=[];
    aux(i,:)=0;
    aux(i,[lst lst+1 lst+2])=1;
    onsets{i}=lst;
end
    
[r,p]=corrcoef(aux');

for i=1:size(aux,1)
    lst=find(p(i,i+1:end)<0.001);
    aux(i+lst,:)=0;
end
[i,j]=find(aux==1);
i=unique(i);

for j=1:length(i)
    if(length(onsets{i(j)})>30)
        st=nirs.design.StimulusEvents;
        st.name=['aux_' num2str(i(j))];
        k=dsearchn(t,onsets{i(j)}'/fs);
        st.onset=t(k);
        st.dur=ones(size(st.onset))*2*mean(diff(t));
        st.amp=ones(size(st.onset));
        stim(st.name)=st;
    end
end

end