function jobs = single_subject
jobs=nirs.modules.ImportData();
jobs.Input='raw';
jobs=nirs.modules.RemoveStimless(jobs);

jobs = nirs.modules.Resample(jobs);
jobs.Fs = 5; % resample to 5 Hz

jobs = nirs.modules.OpticalDensity( jobs );
jobs = nirs.modules.BeerLambertLaw( jobs );
jobs = nirs.modules.ExportData(jobs);
jobs.Output='Hb';
jobs = nirs.modules.TrimBaseline( jobs );
jobs.preBaseline   = 30;
jobs.postBaseline  = 30;
jobs = nirs.modules.AR_IRLS(jobs );
jobs = nirs.modules.ExportData(jobs);
jobs.Output='SubjStats';