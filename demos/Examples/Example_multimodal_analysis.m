% This is an example of mulitmodal (NIRS-EEG) analysis with the toolbox.
% This is mostly work-in-progress

% This simulation function will generate NIRS and EEG data (stored as
% elements in the cell array "data" based on the forward model projection
% of the "truthImage" image.  The truth (cell array) is the expected values
% in channel space

%% 

[data, truth, truthImage] = eeg.testing.simMultimodalData;

eeg_raw=data{1};
nir_raw=data{2};

% Run the NIRS analysis pipeline
job=nirs.modules.Resample;
job=nirs.modules.OpticalDensity(job);
job=nirs.modules.BeerLambertLaw(job);
job=nirs.modules.GLM(job);
NIRS_Stats=job.run(nir_raw);
NIRS_Stats.probe.defaultdrawfcn='10-20';

% Run the EEG analysis pipeline
job=eeg.modules.AverageERP;
job.prewhiten=true;
job.basis=eeg.design.basis.ERP;
EEG_Stats=job.run(eeg_raw);

% Show the results
EEG_Stats.draw('tstat',[],'p<0.05')
NIRS_Stats.draw('tstat',[],'p<0.05')

% Image recon

% For NIRS image recon, we need optical density, so let's jsut rerun
% without the MBLL step
job=nirs.modules.Resample;
job=nirs.modules.OpticalDensity(job);
job=nirs.modules.AR_IRLS(job);
NIRS_Stats=job.run(nir_raw);
NIRS_Stats.probe.defaultdrawfcn='10-20';



Colin=nirs.registration.Colin27.BEM;
Colin.mesh(4)=[];

% EEG forward model
fwdeeg=eeg.forward.FieldTrip;
fwdeeg.mesh=Colin.mesh;
fwdeeg.probe=EEG_Stats.probe;
fwdeeg.prop=[1 NaN NaN];
JacobEEG=fwdeeg.jacobian;

% NIRS forward model
fwdnirs=nirs.forward.ApproxSlab;
fwdnirs.Fm=0;
fwdnirs.mesh=Colin.mesh(3);
lambda=unique(NIRS_Stats.probe.link.type);
fwdnirs.prop=nirs.media.tissues.brain(lambda,.7,50);
fwdnirs.probe=NIRS_Stats.probe.swap_reg;
JacobNIRS=fwdnirs.jacobian('spectral');

% This is the multimodal image recon code, but it can also do
% single-modality too
job=eeg.modules.MultimodalImageReconMFX;
job.mesh=fwdnirs.mesh;
job.mesh.transparency=1;
job.basis = nirs.inverse.basis.gaussian(job.mesh,20);
job.formula = 'beta ~ -1 + cond';  

% EEG
job.jacobian=Dictionary();
job.jacobian('default:eeg')=JacobEEG;
%job.jacobian('default:nirs')=JacobNIRS;
job.probe=Dictionary();
job.probe('default:eeg')=fwdeeg.probe;
%job.probe('default:nirs')=fwdnirs.probe;

prior=[];
prior.eeg=zeros(size(JacobEEG.eeg,2),1);
%prior.hbo=zeros(size(JacobNIRS.hbo,2),1);
%prior.hbr=zeros(size(JacobNIRS.hbr,2),1);
job.prior=Dictionary();
job.prior('default')=prior;


ImageStatsEEG=job.run({EEG_Stats});


% NIRS only
job.jacobian=Dictionary();
%job.jacobian('default:eeg')=JacobEEG;
job.jacobian('default:nirs')=JacobNIRS;
job.probe=Dictionary();
%job.probe('default:eeg')=fwdeeg.probe;
job.probe('default:nirs')=fwdnirs.probe;

prior=[];
%prior.eeg=zeros(size(JacobEEG.eeg,2),1);
prior.hbo=zeros(size(JacobNIRS.hbo,2),1);
prior.hbr=zeros(size(JacobNIRS.hbr,2),1);
job.prior=Dictionary();
job.prior('default')=prior;

ImageStatsNIRS=job.run({NIRS_Stats});

% Multimodal
job.jacobian=Dictionary();
job.jacobian('default:eeg')=JacobEEG;
job.jacobian('default:nirs')=JacobNIRS;
job.probe=Dictionary();
job.probe('default:eeg')=fwdeeg.probe;
job.probe('default:nirs')=fwdnirs.probe;

prior=[];
prior.eeg=zeros(size(JacobEEG.eeg,2),1);
prior.hbo=zeros(size(JacobNIRS.hbo,2),1);
prior.hbr=zeros(size(JacobNIRS.hbr,2),1);
job.prior=Dictionary();
job.prior('default')=prior;

% Multimodal
ImageStatsMM=job.run({NIRS_Stats EEG_Stats});

%NIRS only
ImageStatsNIRS=job.run({NIRS_Stats});

% EEG only
ImageStatsEEG=job.run({EEG_Stats});
