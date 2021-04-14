function data = load_NIRX_processed(file,probe)

load(file);
data = nirs.core.Data;
data.time=[1:size(nirs_data.oxyData,1)]'*nirs_data.samplerate;

if(isstr(probe))
    data.probe=nirs.io.loadNIRxProbe(probe);
else
    data.probe=probe;
end

link=data.probe.link;

data.probe.link=[link; link; link];
data.probe.link.type=[repmat({'hbo'},height(link),1);...
    repmat({'hbr'},height(link),1); repmat({'hbt'},height(link),1)];

nchan=height(link);
data.data=[nirs_data.oxyData(:,1:nchan) nirs_data.dxyData(:,1:nchan) nirs_data.tHbData(:,1:nchan)];

st=nirs.design.StimulusEvents;
st.name='task';
st.onset=nirs_data.onset_task*nirs_data.samplerate;
st.dur=nirs_data.dur_task*nirs_data.samplerate;
st.amp=ones(size(st.onset));

data.stimulus('task')=st;
data.description=file;