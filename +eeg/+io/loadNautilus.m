function data = loadNautilus(filename,fs)
% This function loads g.Nautilus HDF5 EEG data

if(nargin<2)
    fs=250;
end

data=eeg.io.loadHDF5(filename);

fsOrig=data.RawData.DAQDeviceCapabilities.DAQDeviceCapabilities.AnalogChannelProperties.ChannelProperties(1).SampleRate;
d=data.RawData.Samples;

for i=1:length(data.RawData.DAQDeviceCapabilities.DAQDeviceCapabilities.AnalogChannelProperties.ChannelProperties)
    labels{i,1}=data.RawData.DAQDeviceCapabilities.DAQDeviceCapabilities.AnalogChannelProperties.ChannelProperties(i).ChannelName{1};
end

auxinfo=data.AsynchronData.AsynchronSignalTypes.ArrayOfAsynchronSignalDescription.AsynchronSignalDescription;
aux=data.AsynchronData;

tbl=nirs.util.list_1020pts('?');
ismember(lower(tbl.Name),lower(labels));
lst=find(ismember(lower(labels),lower(tbl.Name)));

data=eeg.core.Data;
data.probe=eeg.core.Probe({labels{lst}});
n=round(fsOrig/fs);

d=resample(double(d'),1,n)';
fs=fsOrig/n;


data.data=d(lst,:)';
data.time=[0:size(d,2)-1]/fs;
data.description=which(filename);

stim=Dictionary;
ut=unique(aux.TypeID);
for i=1:length(ut)
    st=nirs.design.StimulusEvents;
    ID=find(vertcat(auxinfo.ID)==ut(i));
    st.name=auxinfo(ID).Name;
    st.onset=aux.Time(find(aux.TypeID==ut(i)))/fsOrig;
    st.amp=aux.Value(find(aux.TypeID==ut(i)));
    st.dur=ones(size(st.onset))*2/fs;
    stim(st.name)=st;
end


data.stimulus=stim;
data.description=filename;
end
