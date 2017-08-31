function data = add_MRI_pulsetiming(data,chanIdx)

if(nargin<2)
    chanIdx=2;
end
for i=1:length(data)
    [aux,t]=eeg.io.loadBrainVisionAux(data(i).description);
    lst=find(aux(chanIdx,:)>max(aux(chanIdx,:))*.8);
    lst(1+find(diff(lst)<10))=[];
    st=nirs.design.StimulusEvents;
    st.onset=t(lst)';
    st.name='MRI pulse';
    st.dur=0.1*ones(size(st.onset));
    st.amp=ones(size(st.onset));
    data(i).stimulus('MRI_Pulse')=st;
end

return
