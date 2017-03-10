function stim = loadFiffStimulus(filename)
% This function will read the stimulus info from MEG Fiff files

if(iscell(filename))
    for i=1:length(filename)
        stim(i)=eeg.io.loadFiffStimulus(filename{i});
    end
    return
end


%hdr=fiff_read_meas_info(filename);
raw=fiff_setup_read_raw(filename);
hdr=raw.info;
aux=fiff_read_raw_segment(raw,0,1E6,find(cat(hdr.chs.kind)==3));
time=[0:size(aux,2)-1]'/hdr.sfreq;
stim=findstim(aux,hdr.sfreq,time);

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