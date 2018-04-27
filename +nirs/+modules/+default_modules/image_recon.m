function jobs = image_recon
% This pipeline will run the simple image reconstruction 
% the raw data must have a Probe1020 type object and is tested in 
% the asset step (job step #2)



jobs=nirs.modules.ImportData();
jobs.Input='raw';

jobs=nirs.modules.Assert(jobs);
jobs.throwerror=true;
jobs.condition=@(data)isa(data.probe,'nirs.core.Probe1020');

jobs=nirs.modules.RemoveStimless(jobs);
jobs = nirs.modules.FixNaNs(jobs);
jobs = nirs.modules.Resample(jobs);
jobs.Fs = 5; % resample to 5 Hz

jobs = nirs.modules.OpticalDensity( jobs );
jobs = nirs.modules.ExportData(jobs);
jobs.Output='dOD';
jobs = nirs.modules.TrimBaseline( jobs );
jobs.preBaseline   = 30;
jobs.postBaseline  = 30;
jobs = nirs.modules.AR_IRLS(jobs );
jobs = nirs.modules.ExportData(jobs);
jobs.Output='SubjStats';

jobs = nirs.modules.ImageReconMFX(jobs);

Slab = nirs.forward.ApproxSlab;
probe=evalin('base','raw(1).probe');
lambda=unique(probe.link.type);
Slab.prop=nirs.media.tissues.brain(.7,50,lambda);
Slab.mesh=probe.getmesh();
Slab.mesh=Slab.mesh(end);
Jac=Slab.jacobian('spectral');
jobs.probe('default')=Slab.probe;
jobs.jacobian('default')=Jac;
jobs.formula='beta ~ -1 + cond';
jobs.mesh=Slab.mesh;

jobs.basis=nirs.inverse.basis.identity(size(Slab.mesh.nodes,1));

%jobs.basis=nirs.inverse.basis.gaussian(Slab.mesh,30);

jobs = nirs.modules.ExportData(jobs);
jobs.Output='ImageStats';



