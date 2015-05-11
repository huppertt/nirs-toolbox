
data = nirs.io.loadDirectory( '~/resting_state', {} );

% covert to optical density
j = nirs.modules.OpticalDensity( );
j = nirs.modules.BeerLambertLaw( j );

% resample to 4 Hz
j = nirs.modules.Resample( j );
j.Fs = 4;

% subject level
j = nirs.modules.AR_IRLS( j );

% ROC
r = nirs.testing.SubjectLevelROC( j );
r.beta  = 4; % uM
r.niter = 200;

r = r.run( data );