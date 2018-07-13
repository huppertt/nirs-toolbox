% Script to run EEG analysis of Balance data

EEGfiles{1}='test_polhemus5.vhdr';
EEGfiles{2}='test_polhemus6.vhdr';
EEGfiles{3}='test_polhemus7.vhdr';
EEGfiles{4}='test_polhemus8.vhdr';

EEGtasks{1}='SOT1->SOT4';
EEGtasks{2}='SOT1->SOT4';
EEGtasks{3}='SOT2->SOT5';
EEGtasks{4}='SOT2->SOT5';

% I don't have a good notation for stim marks in the EEG files yet
for i=1:length(EEGfiles)
    raw(i) = eeg.io.loadBrainVision(EEGfiles{i},500);
    endtime = raw(i).time(end);
    stim=nirs.design.StimulusEvents;
    stim.name=EEGtasks{i};
    stim.onset=endtime-120;
    stim.dur=60;
    stim.amp=1;
    raw(i).stimulus(EEGtasks{i})=stim;
end

%% Run the preprocessing
job=eeg.modules.BandPassFilter;
job.lowpass=[];
job.highpass=1;
job=eeg.modules.KurtoisFilter(job);
job=nirs.modules.Resample(job);
job.Fs=100;
job=nirs.modules.PCAFilter(job);
job.ncomp=.4;

job=eeg.modules.BandPassFilter(job);
job.lowpass=25;
job.highpass=1;

processed=job.run(raw);


%% Frequency analysis
job=eeg.modules.WaveletTransform;
job.scale_smoothing=true;
job.frequencies(end,:)=[];  % get rid of Gamma (25-50hz) 
FreqData = job.run(processed);

%% Now, run the block average analysis model
job=eeg.modules.AverageERP();
job.prewhiten=true;
job.basis=nirs.design.basis.Canonical;
FDStats=job.run(FreqData);


% Group level model
job=eeg.modules.MixedEffects;
job.formula = 'beta ~ -1 + cond:freq'; 

GroupStats=job.run(FDStats);


% Image recon

% Do an image reconstruction of the results
job=eeg.modules.ImageReconMFX;

FTfwd=eeg.forward.FieldTrip;
FTfwd.mesh=nirs.registration.Colin27.BEM;
FTfwd.probe=FDStats(1).probe;
FTfwd.prop=[1 NaN 1 1];

Jacob=FTfwd.jacobian;
job.jacobian('default')=Jacob;   
job.probe('default')=processed.probe;
job.mesh=FTfwd.mesh.mesh(end);  

job.formula = 'beta ~ -1 + cond:freq';  % Simple fixed effects model
%job.basis=nirs.inverse.basis.identity(job.mesh);
job.basis = nirs.inverse.basis.gaussian(job.mesh,20);

% Now create the priors in the model
job.mask =[];
% This is the Minimum Norm estimate
prior.eeg=zeros(size(Jacob.eeg,2),1);
job.prior=Dictionary();
job.prior('default')=prior;


ImageStats=job.run(FDStats);

% To draw
ImageStats.draw('tstat',[],'p<0.05',[],{'left','right'})