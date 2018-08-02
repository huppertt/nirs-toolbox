function data = loadFiff(filename,fs)
% This function will read MEG Fiff files

if(nargin<2)
    fs=250;
end

if(iscell(filename))
    for i=1:length(filename)
        data(i)=eeg.io.loadFiff(filename{i},fs);
    end
    return
end


%hdr=fiff_read_meas_info(filename);
raw=fiff_setup_read_raw(filename);
hdr=raw.info;
d=fiff_read_raw_segment(raw);
a={hdr.chs.kind};
for i=1:length(a); kind(i)=a{i}; end;
aux=d(find(kind==3),:);

if(fs<hdr.sfreq/2)
    n=round(hdr.sfreq/fs);
    d=resample(double(d'),1,n)';
    fs=hdr.sfreq/n;
else
    fs=hdr.sfreq;
end


data=eeg.core.Data;
lst=find([hdr.chs(:).kind]==1);
data.data=d(lst,:)';
data.time=[0:size(d,2)-1]/fs;
data.description=which(filename);

data.probe=eeg.core.MEGProbe(hdr);
try
stim=findstim(aux,hdr.sfreq,data.time);
data.stimulus=stim;
end
data.description=filename;
end

function stim = findstim(aux,fs,t)

aux2=aux;
stim=Dictionary;
for i=1:size(aux,1)
    s=diff(aux(i,:));
    s=s-s(1);
    s=s./sqrt(var(s));
    lst=find(s>5);
    lst(find(diff(lst)<50))=[];
    aux(i,:)=0;
    aux(i,[lst lst+1 lst+2])=1;
    onsets{i}=lst;
    
    lst=find(s<-5);
    if(~isempty(lst))
        lst(find(diff(lst)<50))=[];
        offsets{i}=lst;
    else
        offsets{i}=[];
    end
end
offsets{end+1}=t(end);

[r,p]=corrcoef(aux');

for i=1:size(aux,1)
    lst=find(p(i,i+1:end)<0.001);
    aux(i+lst,:)=0;
end
[i,j]=find(aux==1);
i=unique(i);

for j=1:length(i)
    if(length(onsets{i(j)})>0)
        st=nirs.design.StimulusEvents;
        st.name=['aux_' num2str(i(j))];
        k=dsearchn(t,onsets{i(j)}'/fs);
        k1=dsearchn(t,offsets{i(j)}'/fs);
        st.onset=t(k);
            tempmin=min(length(t(k)), length(t(k1)));
            k(tempmin+1:end)=[];
            k1(tempmin+1:end)=[];
        st.dur=t(k1)-t(k);
        st.amp=ones(size(st.onset));
        stim(st.name)=st;
    end
end

end