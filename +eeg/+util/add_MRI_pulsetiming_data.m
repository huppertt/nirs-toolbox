function data = add_MRI_pulsetiming_data(data,TR)

for i=1:length(data)
    
    [hdr.fs, hdr.label,hdr.meta]=bva_readheader(data(i).description);
    d=bva_loadeeg(data(i).description);
    t=[0:size(d,2)-1]./hdr.fs;
    
    d=mean(d,1);
    
    lst=find(d>max(d)*.8);
    lst(1+find(diff(lst)<2))=[];
    lstkeep=1;
    
    trs=mode(diff(t(lst)));
    TR=round(TR/trs)*trs;
    
    while(1)
        tt=t(lst(lstkeep(end)));
        ii=min(find(t(lst)>=(tt+TR)));
        if(isempty(ii))
            break
        end
        lstkeep(end+1)=ii;
    end
    lst=lst(lstkeep);
    
    
    st=nirs.design.StimulusEvents;
    st.onset=t(lst)';
    st.name='MRI pulse';
    st.dur=0.1*ones(size(st.onset));
    st.amp=ones(size(st.onset));
    data(i).stimulus('MRI_Pulse')=st;
    
end