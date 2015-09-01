%% This is a demo to show how to import *.nirs data from the HOMER-2 format 
% The HOMER-2 format uses the field "CondNames" and the "s" variable to
% define stim timing, which differs from the StimDesign or stimulus
% Dictionary formats used in this tool box.  Namely, the HOMER-2 format does not 
% specify duration or amplitude information which is needed for the GLM
% model and so these needs to be added.

% Data (in the HOMER-2) format is loaded in the same way
raw = nirs.io.loadDirectory(folder,{'group','subject'});

% The load function will use the 's' variable and the CondNames to define
% the stimulus info.  
% If the "s" variable is blocked (e.g. s channel is 1 for the entire
% duration of the stim rather then just the onset), then the duration will
% be fine and the next step can be skipped.  Otherwise (e.g. only the onsets are denoted),
% we have to change the duration from a delta function (single point) to the 
% duration of the task.

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

% To chnage the timing of stimulus events you modify the entry of the
% stimTable using NaN as placeholds to indicate keeping the old values

% E.g to change the duration of the events to 10s but keep the onsets as
% defined by the data import function 

% This will cut off all but the first entry of the table
stimTable(2:end,:)=[];
stimTable.FileIdx=NaN;  % If we set the FileIdx to NaN then it will apply to all files.
% This is a short-cut which does the same as leaving the entire table
% intact and changing every entry seperately.  If we left this as
% FileIdx=1, then our change would only be applied to the first file

stimTable.Talk.onset=NaN;  % Using NaN will keep the original timing.  Again, if we left the 
% values as real values (which in this case correspond to the first file),
% then it will replace all the times (and potentially mess up files 2 and
% beyond).  
stimTable.Talk.dur=10; % But let's change the duration to 10s for all trials

% If there are conditions in the stimTble that we don't want to change at
% all, then we can just remove them from the table and no changes will be
% done to these variables in the data file (e.g. same as using NaN for all
% fields in the stimulusEvent variable).
stimTable.Rest=[];

% Likewise, we can add new conditions by adding new columns to the table
% Let's try a stimusVector object

stimTable.NewCondition=nirs.design.StimulusVector;
stimTable.NewCondition.name='NewCondition';
stimTable.NewCondition.time=[0:1:300];
stimTable.NewCondition.vector = sin(stimTable.NewCondition.time);


% To implement the changes, we use the ChangeStimulusInfo job
j = nirs.modules.ChangeStimulusInfo();
j.ChangeTable = stimTable;

rawChanged = j.run(raw);

% if we draw it, we can see the changes
% either by
rawChanged(1).draw;
% or using the nirsviewer GUI

nirs.viz.nirsviewer;



