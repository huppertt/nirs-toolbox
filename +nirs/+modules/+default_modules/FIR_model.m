function jobs = FIR_model

Fs=1;

jobs=nirs.modules.ImportData();
jobs.Input='raw';
jobs=nirs.modules.RemoveStimless(jobs);
jobs = nirs.modules.Resample(jobs);
jobs.Fs = Fs; % resample to 1 Hz
 jobs = nirs.modules.OpticalDensity( jobs );
jobs = nirs.modules.BeerLambertLaw( jobs );
jobs = nirs.modules.ExportData(jobs);
jobs.Output='Hb';
jobs = nirs.modules.TrimBaseline( jobs );
jobs.preBaseline   = 30;
jobs.postBaseline  = 30;
jobs = nirs.modules.AR_IRLS(jobs );
jobs.goforit=true;
bas=nirs.design.basis.FIR;
bas.isIRF=true;
bas.nbins=16*Fs;
bas.binwidth=1;
jobs.basis('default')=bas;

jobs = nirs.modules.ExportData(jobs);
jobs.Output='SubjStats';
jobs = nirs.modules.GroupAverage(jobs );
jobs.formula       = 'beta ~ -1 + cond';  % See help fitlme for examples
jobs = nirs.modules.ExportData(jobs);
jobs.Output='GroupStats';

