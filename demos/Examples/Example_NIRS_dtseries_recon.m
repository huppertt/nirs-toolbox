% Example of creating "brain-space" NIRS time-course data

raw =nirs.testing.simData_registered;

job=nirs.modules.OpticalDensity;
job=dtseries.modules.Convert2dtseries(job);

HbBrain = job.run(raw);

% 
% HbBrain = 
% 
%   Data with properties:
% 
%      description: []
%             data: [3001×32 double]
%             mesh: [1×1 dtseries.core.Mesh]
%             time: [3001×1 double]
%       projectors: [9698×16 double]
%              cov: [32×32 double]
%         stimulus: [1×1 Dictionary]
%     demographics: [1×1 Dictionary]
%               Fs: 10

% Note- 
%  1) data is stored as its eigenvectors and must be "projected" to 
%     get into the actual brain space 
%   you can use:
%       data = HbBrain.getdata;  
%   to do this.  This is true of the covariance model too
%  2) Hbbrain.mesh is roughly equivelent to the "probe" field in the other
%     classes

