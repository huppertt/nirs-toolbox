%% This is a demo to show how to import *.nirs data from the HOMER-2 format 
% The HOMER-2 format uses the field "CondNames" and the "s" variable to
% define stim timing, which differs from the StimDesign or stimulus
% Dictionary formats used in this tool box.  Namely, the HOMER-2 format does not 
% specify duration or amplitude information which is needed for the GLM
% model and so these needs to be added.

clear 
% change this to save results somewhere else
root_dir = '/Users/thuppert/Desktop/tmp' ;


% Download the sample data from the HOMER-2 site
if(~exist(root_dir,'dir') || ~exist(fullfile(root_dir,'SampleData'),'dir'))
    mkdir(root_dir);
    disp('downloading sample data from HOMER-2 website');
    %% download the dataset
    urlwrite('http://www.nmr.mgh.harvard.edu/optics/resources/homer2/SampleData.zip', ...
        [root_dir filesep 'demo_homer2.zip']);
    % This command will download the demo_data.zip file from the server.  This
    % step can be skipped if you already downloaded this. This could take a few minutes if your internet conenction is slow
    % The file is about 90Mb in size.
    
    % unzip the data
    unzip([root_dir filesep 'demo_homer2.zip'],[root_dir filesep]);
       
else
    disp(['Data found in: ' root_dir ': skipping download']);
end


demo = 'Example1_Simple_Probe';

% Data (in the HOMER-2) format is loaded in the same way
raw = nirs.io.loadDirectory(fullfile(root_dir,'SampleData',demo),{});
% The load function will use the 's' variable and the CondNames to define
% the stimulus info.  
% If the "s" variable is blocked (e.g. s channel is 1 for the entire
% duration of the stim rather then just the onset), then the duration will
% be fine and the next step can be skipped.  Otherwise (e.g. only the onsets are denoted),
% we have to change the duration from a delta function (single point) to the 
% duration of the task.

% In this file, the stim marks were not preprocessed, so the code read the
% timing from the "s" field and used a default name (stim_condition1)

% This function will extract all the stim info from the data variable
stimTable=nirs.createStimulusTable(raw);
disp(stimTable);

% The height of this table matches the length of the "raw" data field.
% The first column denotes the FileIdx.
% The FileIdx is used to match the stiminfo to the data.  This allows us to
% modifiy each file individually, only certain files (if you remove all but 
% these rows from the table), or modifiy all the stim conditions together.
% We use "NaN" as a wildcard varible.  Thus, if FileIdx=NaN, this entry
% will be applied to all the files.  

% The entries in the table for each condition are either type
% nirs.design.StimulusEvents or nirs.design.StimVector;
%
% nirs.design.StimulusEvents - denote indidual events and have the fields
%      name: - the name of the condition
%     onset: - a vector of onset times (in seconds)
%       dur: - a vector of durations for each event
%       amp: - the amplitude of each event (used in parametric designs)
%
% In this notation, the duration for eact trial can be independently set.
% If a scalar (single value) is given for either amp or dur, then it is
% assumed that all trials have the same value (and differ only in onset time).
%
% nirs.design.StimulusVector- denotes a continious regressor and has the fields
%       name: - the name of the condition
%     vector: - a vector of amplitudes (per time)
%       time: - the time vector for the samples.
%
%  Note- The spline-interp function is used to resample the vector onto the data
%  when creating the design matrix

% To chabge the timing of stimulus events you modify the entry of the
% stimTable using NaN as placeholds to indicate keeping the old values

% E.g to change the duration of the events to 10s but keep the onsets as
% defined by the data import function 

% This will cut off all but the first entry of the table
stimTable(2:end,:)=[];  
stimTable.FileIdx=NaN;  % If we set the FileIdx to NaN then it will apply to all files.
% This is a short-cut which does the same as leaving the entire table
% intact and changing every entry seperately.  If we left this as
% FileIdx=1, then our change would only be applied to the first file

stimTable.stim_channel1.onset=NaN;  % Using NaN will keep the original timing.  Again, if we left the 
% values as real values (which in this case correspond to the first file),
% then it will replace all the times (and potentially mess up files 2 and
% beyond).  
stimTable.stim_channel1.dur=5; % But let's change the duration to 10s for all trials


% To implement the changes, we use the ChangeStimulusInfo job
j = nirs.modules.ChangeStimulusInfo();
j.ChangeTable = stimTable;

% Now, run the actual job
rawChanged = j.run(raw);


% if we draw it, we can see the changes
% either by
rawChanged(1).draw;
% or using the nirsviewer GUI
% >> nirs.viz.nirsviewer;

% Now, let's run an analysis on this

% Create the pipeline (this is the most basic pipeline possible)
j = nirs.modules.OpticalDensity();
j = nirs.modules.BeerLambertLaw(j);
j = nirs.modules.AR_IRLS(j);

% Let's use a FIR model to show how to do deconvolution
FIRbasis=nirs.design.basis.FIR;

% FIRbasis = 
%   FIR with properties:
%     nbins: 10  - number of bins in the HRF 
%     binwidth: 5 - number of time-points per bin

% Currently our sample rate is 10Hz (=raw(1).Fs)
% Thus, binwidth =5 corresponds to a 0.5s bin
% The nbins =10 tells us the window will be 0.5s x 10 = 5s wide, which is not
% wide enough.  Note, this is for the estimation of the impulse response
% function (so the response is convolved with whatever the stimulus
% duration is).  Let's make the width 15s by setting nbins =30
FIRbasis.nbins=30;

j.basis('default')=FIRbasis;

SubjStats = j.run(rawChanged);
% Note, SubjStats now has 240 betas (8 channels x 30 betas per condition)
% If we used SubjStats.draw() it would draw 30 images (one per time point in the HRF)

% This funciton will create the HRF data from the basis set.  This function
% also needs to know the duration of the event (default = impulse response)
% and the sample rate.
HRF = nirs.design.extractHRF(SubjStats,j.basis,5,rawChanged.Fs);

% This is a core.Data object so we can draw it (or load it in
% nirs.viz.nirsviewer
HRF.draw();

% To make a contast over some time window
c = zeros(1,30);
c(1,[8:16])=1;  % Set the contrast from 4s-8s (remember this the impulse response not the HRF per se). 
Contrast = SubjStats.ttest(c);

Contrast.draw();





