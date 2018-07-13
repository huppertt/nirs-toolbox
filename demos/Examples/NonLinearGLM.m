

% THIS DOES NOT WORK...   Its in progress.

% data = nirs.testing.simData;
% 
% job=nirs.modules.Resample();
% job.Fs=1;
% job=nirs.modules.OpticalDensity(job);
% job=nirs.modules.BeerLambertLaw(job);
% job=nirs.modules.nonlin_GLM(job);

num_iter=5;

ROCtest=nirs.testing.ChannelStatsROC;
ROCtest.simfunc=@()nirs.testing.simData;

jobs = nirs.modules.default_modules.single_subject;
list=nirs.modules.pipelineToList(jobs);
ROCtest.pipeline=nirs.modules.listToPipeline(list);

list{8}=nirs.modules.nonlin_GLM;
ROCtest.pipeline(2)=nirs.modules.listToPipeline(list);


ROCtest=ROCtest.run(num_iter);
ROCtest.draw;



num_iter=50;

ROCtest=nirs.testing.ChannelStatsROC;
ROCtest.simfunc=@()nirs.testing.simDataWithSuperficial;

jobs = nirs.modules.default_modules.single_subject;
list=nirs.modules.pipelineToList(jobs);
ROCtest.pipeline=nirs.modules.listToPipeline(list);

list{end+1}=nirs.modules.TestSuperficial;
ROCtest.pipeline(2)=nirs.modules.listToPipeline(list);
list{end}.method=1;
ROCtest.pipeline(3)=nirs.modules.listToPipeline(list);

ROCtest=ROCtest.run(num_iter);
ROCtest.draw;

