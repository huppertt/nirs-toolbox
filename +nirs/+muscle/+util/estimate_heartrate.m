function estimate_heartrate(hb)

% first filter the date
j=eeg.modules.BandPassFilter;

j.highpass=.1;
j.lowpass=min([5 hb(1).Fs/2-.05*hb.Fs/2]);

hbf=j.run(hb);

hbo=mean(hbf.data(:,ismember(hbf.probe.link.type,'hbo')),2);
hbr=mean(hbf.data(:,ismember(hbf.probe.link.type,'hbr')),2);