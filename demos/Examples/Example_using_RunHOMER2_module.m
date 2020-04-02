
% create a test scan with 20s on/off block events
noise = nirs.testing.simARNoise;
stim = nirs.testing.blockedStimDesign(noise.time,20,20);
raw=nirs.testing.simData(noise,stim);



job=nirs.modules.OpticalDensity;
job=nirs.modules.BeerLambertLaw(job);
job=nirs.modules.Run_HOMER2(job);
job.fcn='hmrBlockAvg';   

% note- this module populates all the input/outputs automatically.  Any
% HOMER2 variable (e.g. nTrials for the hmrBlockAvg function) gets plaes in
% a tmp field which can be accessed by any other HOMER function in the
% pipeline.  Note the structure of the input and output variable handling
% in this function: 
% job.outputs
%   6×2 cell array
%     {'yavg'   }    {@(data,d)setfield(data,'data',d)}
%     {'ystd'   }    {                 @(data,tmp)data}
%     {'tHRF'   }    {@(data,t)setfield(data,'time',t)}
%     {'nTrials'}    {                 @(data,tmp)data}
%     {'ysum2'  }    {                 @(data,tmp)data}
%     {'yTrials'}    {                 @(data,tmp)data}
% 
% This controls the interfacing between my code and HOMER2

% since this blcok avg function has a  "trange" input that needs user
% input, you need to edit the job.vars  fields (which are automatically
% populated if possible)
job.vars.trange=[-5 20];

HRF = job.run(raw);