clear 

root_dir = 'set this' ;

%% load data
raw = nirs.io.loadDirectory(root_dir, {'group', 'subject'});

%% preprocessing pipeline

% remove junk files 
j = nirs.modules.RemoveStimless( );

% resample to 4 Hz
j = nirs.modules.Resample( j );
j.Fs = 4;

% limit leading and post baseline to 60s
j = nirs.modules.TrimBaseline( j );
j.preBaseline  = 60;
j.postBaseline = 60;

% covert to optical density
j = nirs.modules.OpticalDensity( j );

% convert to hemoglobin
j = nirs.modules.BeerLambertLaw( j );

hb = j.run( raw );


%% subject level pipeline
j = nirs.modules.AR_IRLS();
j.verbose = true;

j.trend_func = @(t) nirs.design.trend.legendre(t, 3);

S = j.run( hb );

%% group level
clear j
j = nirs.modules.MixedEffects( );
j.formula = 'beta ~ -1 + group:cond + (1|subject)';
j.dummyCoding = 'full';

G = j.run( S );