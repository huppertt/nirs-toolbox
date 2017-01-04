num_iter=5;

% THIS DOES NOT WORK...   Its in progress.


nlbas = nirs.design.basis.nonlinearHRF();

%Simulate data
nirs_location = '/Volumes/sparky/R21/R21 DATA/Brainy Kids Data for Ted/Processed Brainy Kids Data/';
raw = nirs.io.loadDirectory(nirs_location,{'subject'});



job=nirs.modules.RemoveShortScans;
raw=job.run(raw);

can = nirs.design.basis.Canonical;
%TODO - change from the defaults;

basis=Dictionary;
basis('default')=can;
BetaAmp=7;
ROCtest=nirs.testing.ChannelStatsROC;
ROCtest.simfunc=@()nirs.testing.simData(raw(randi(length(raw),1,1)),...
    [],BetaAmp,[],basis);

jobs = nirs.modules.default_modules.single_subject;
list=nirs.modules.pipelineToList(jobs);
ROCtest.pipeline=nirs.modules.listToPipeline(list);

list{8}=nirs.modules.nonlin_GLM;
ROCtest.pipeline(2)=nirs.modules.listToPipeline(list);





ROCtest=ROCtest.run(num_iter);
ROCtest.draw;