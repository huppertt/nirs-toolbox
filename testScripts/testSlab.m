clear

% segmented image
nirsImage = nirs.Image( uint8(ones(200,200,100)) );

% make a probe
detPos = [100 100 1];
srcPos = [110 100 1; 120 100 1; 130 100 1; 110 100 1; 120 100 1; 130 100 1];
srcDir = repmat([0 0 1],[size(srcPos,1) 1]);
detDir = repmat([0 0 1],[size(detPos,1) 1]);

link = [1 1 1;
    2 1 1;
    3 1 1;
    4 1 2;
    5 1 2;
    6 1 2];

lambda = [690; 830];

nirsProbe = nirs.Probe( srcPos, detPos, link, lambda, [], [], srcDir, detDir );

% define optical properties
optProp = nirs.OpticalProperties( [.01 .03]', [1 1.2]', [1.3 1.3]', [690 830]' );

% slab forward model
slabModel = nirs.SlabForwardModel( nirsProbe, optProp );

meas = slabModel.measurement();
estProp = nirs.slabSolver( meas )

% create forward model object
fwdModel = nirs.MCXForwardModel( nirsImage, optProp, nirsProbe );
% fwdModel.modFreq = 110e6;
% fwdModel.nPhotons = 1e7;
% fwdModel.timeStep = 1/110e6/32;
% fwdModel.nTimeGates = 32;

% measurement
meas = fwdModel.measurement();

% test with slab model
estProp = nirs.slabSolver( meas )