
data = nirs.io.loadDirectory( '~/resting_state', {} );

%% AR-IRLS
% covert to optical density
j = nirs.modules.OpticalDensity( );
j = nirs.modules.BeerLambertLaw( j );

% resample to 4 Hz
j = nirs.modules.Resample( j );
j.Fs = 4;

% subject level
j = nirs.modules.AR_IRLS( j );

% ROC
r1 = nirs.testing.SubjectLevelROC( j );
r1.beta  = 5; % uM
r1.niter = 10;

r1 = r1.run( data );

%% Multivariate AR-IRLS
% covert to optical density
j = nirs.modules.OpticalDensity( );

% resample to 4 Hz
j = nirs.modules.Resample( j );
j.Fs = 4;

% subject level
j = nirs.modules.MVGLM( j );

% ROC
r2 = nirs.testing.MVSubjectLevelROC( j );
r2.beta  = 8; % uM
r2.niter = 50;

r2 = r2.run( data );