% Demo of basic EEG analysis with the toolbox
% In this demo, I will show how to load, preprocess, 
% average, and image reconstruct EEG data 
%
% Disclaimer- Currently, I only have code for a brainvision
% EEG system and our focus is multimodal (NIRS-EEG) analysis.
% If you are just looking for EEG (not multimodal), you might be 
% better off with EEGLab or FieldTrip since most of my EEG code
% just wraps their functionality into my toolbox format.

% The +EEG tools are sructured the same way as the base +NIRS toolbox
% including the structure of the data and statistics object classes and
% most of the module and methods functions to interacte with these classes.
% In most cases, the EEG modules and NIRS modules will work on
% both types of data (although obviously running the Beer-Lambert law
% on EEG data is probably not advisable).  But the resample, bandpass
% filtering, PCA, etc code might make sense in either modality

% TODO- In the not so distant future, I will add support for MEG as well in
% this toolbox.

%% Loading data
% Data is loaded via the +eeg/+io functions.  Currently, I only have code
% to load BrainVision EEG files one at a time, but I will soon add more
% data/instrument types and a LoadDirectory function.

% The load function takes a second argument (resample date) which will
% automatically resample the data upon loading.  This is done to save
% memory when loading large data sets.  This is the same as loading and
% then resampling in two job steps

data(1)=eeg.testing.simData;
data(2)=eeg.testing.simData;


% data(1)=eeg.io.loadBrainVision('ted_standing_opticflow.vhdr',200);
% data(2)=eeg.io.loadBrainVision('ted_standing_random_dots.vhdr',200);

% The data class (eeg.core.Data) has the same basic structure as the NIRS
% data class:
%   Data with properties:
%      description: 'ted_standing_opticflow.vhdr'
%             data: [69272x31 double]
%            probe: [1x1 eeg.core.Probe]
%             time: [69272x1 double]
%         stimulus: [1x1 Dictionary]
%     demographics: [1x1 Dictionary]
%               Fs: 200

% This also has the same methods included in the NIRS object including draw
% and sort.  For example, 

data(1).draw(1);  % Will draw the first channel of EEG data

% Stimulus and demographics information are the same as in the NIRS data
% class.

% The eeg.core.Probe class (akin to the nirs.core.Probe class) is used to
% define the identities of the (in this case) 31 EEG channels.
% The EEG probe has few fields then the NIRS equivelent
%
%   Probe with properties:
%     electrodes: [31x5 table]
%           link: [31x2 table]
%
% electrodes is a table (e.g.)
data(1).probe.electrodes
%      Name        X           Y           Z        Units
%     ______    ________    ________    ________    _____
%     'Fp1'      -27.072       100.7      -6.615    'mm' 
%     'Fp2'       32.066      99.815     -10.592    'mm' 
%     'F7'       -70.062      61.802     -18.346    'mm' 
%     'F3'       -47.065       62.85      41.026    'mm' 
%     'Fz'        4.4927      62.878      65.607    'mm' 
%     'F4'        55.162      61.353      33.356    'mm' 

% The 3D positions come from the 10-20 positions on the Colin-27 atlas.
% TODO- add code to resize the head and re-register to landmarks

% The link variable (similar to its use in the NIRS to denote
% source-detectors) is simple at the moment.  This will be modified when we
% add frequency domain analysis later in this demo.
data(1).probe.link
%     electrode    type 
%     _________    _____
%      1           'eeg'
%      2           'eeg'
%      3           'eeg'
%      4           'eeg'
%      5           'eeg'

% The probe can also draw itself.  Because the probe is already registered 
% to the 10-20 system, it will draw on a topo-map by default.  This is akin
% to the registered (nirs.core.Probe1020) data class

data(1).probe.draw


%% Adding stimulus information
% The EEG data re-uses the exact StimulusEvents and StimulusVector classes
% used in the NIRS toolbox.  Here, we add to blocks of 60s long events to
% the data.  
stim=nirs.design.StimulusEvents;
stim.onset=[90 210];
stim.dur=[60 60];
stim.amp=[1 1];
stim.name='Vection';
data(1).stimulus('Vection')=stim;

stim=nirs.design.StimulusEvents;
stim.onset=[60 150];
stim.dur=[60 60];
stim.amp=[1 1];
stim.name='Random';
data(2).stimulus('Random')=stim;

% Now the stimulus information shows in the draw function
data(1).draw(1)

% Demographics information can be added in the same way 
data(1).demographics('subject')='S1';
data(1).demographics('age')=3;
data(2).demographics('subject')='S2';
data(2).demographics('age')=5;

% and is shown usuing the function 
nirs.createDemographicsTable(data)

%     subject    age
%     _______    ___
%     'S1'       3  
%     'S2'       5  

% The nirs.modules.AddDemographics will also work to add demographics from
% a excell sheet (see fnirs_analysis_demo)

%% Basic processing
% This the NIRS toolbox, the EEG tools use a job/module format.  Many of
% the jobs from the NIRS toolbox (e.g. resample and filter) are just reused
% for the EEG (I am not sure if this is lazy or clever).

% For example
j=eeg.modules.KurtoisFilter;  % This is an EEG job to remove eye-blinks and artifacts
                              % This is a PCA or ICA filter followed by
                              % downweighting (via bisquare) of components with
                              % high kurtois compared to the remaining
                              % channels
                              
j=nirs.modules.PCAFilter(j);  % This is reused the NIRS PCA filter for removing spatial covariance

j=eeg.modules.BandPassFilter(j);  % This is a bandpass filter.  Apperently, I 
                                  % had not written a NIRS filter yet.  This module will work fine to
                                  % use on nirs data as well.
j.lowpass=25;   % passbands for the filter
j.highpass=1;   % note, the bandpass filter also resamples to the Nyquist if a lowpass is used


% As with the NIRS data, you create jobs or processing streams by chaining
% jobs together and then issue the run command to process the data
data=j.run(data);

%% Frequency domain analysis
% This particular data is from a NIRS-EEG study on balance using a 60s
% task block.  It doesn't make sense to use a conventional ERP style analysis and instead let's
% look at frequency-domain changes in alpha, beta, etc bands

% The Wavelet transform module converts to frequency data
j=eeg.modules.WaveletTransform;
% This uses the Matlab CWT function and requires the Matlab Wavelet Toolbox 

%   WaveletTransform with properties:
%               wname: 'morl'         <-- wavelet type (anything the cwt
%                                         functon supports)
%     scale_smoothing: 1              <-- preform smoothing
%         frequencies: [6x3 table]
%       convert2power: 1              <-- returns power vs amp/phase data
%                name: 'Wavelet Transform using CWT'

% The frequecies table defines which bands are recovered and can be
% modified.  This is the default values used (in Hz)
% j.frequencies
%     lower    upper     name  
%     _____    _____    _______
%        0        4     'delta'
%        6        7     'theta'
%      7.5     12.5     'alpha'
%        9       11     'beta' 
%     12.5       30     'mu'   
%       25       50     'gamma'

FD = j.run(data);
% This is not terribly slow (~1min per file) and will give a progress update to the screen.




% Run the block averaging model on the Freq-Domain results
j=eeg.modules.AverageERP();
j.prewhiten=true;
j.basis=nirs.design.basis.Canonical;
FDStats=j.run(FD);



% Do an image reconstruction of the results
j=eeg.modules.ImageReconMFX;

FTfwd=eeg.forward.FieldTrip;
FTfwd.mesh=nirs.registration.Colin27.BEM;
FTfwd.probe=FDStats(1).probe;
FTfwd.prop=[1 NaN 1 1];

J=FTfwd.jacobian;
j.jacobian('default')=J;   
j.probe('default')=data.probe;
j.mesh=FTfwd.mesh.mesh(end);  

j.formula = 'beta ~ -1 + cond:freq';  % Simple fixed effects model
%j.basis=nirs.inverse.basis.identity(j.mesh);
%j.basis = nirs.inverse.basis.gaussian(j.mesh,20);


ImageStats=j.run(FDStats);



