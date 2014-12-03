clear

% segmented image
nirsImage = nirs.Image( uint8(ones(200,200,100)) );

% make a probe
detPos = [70 100 1];
srcPos = [80 100 1; 85 100 1; 90 100 1; 95 100 1; 100 100 1];
srcDir = repmat([0 0 1],[size(srcPos,1) 1]);
detDir = repmat([0 0 1],[size(detPos,1) 1]);

link = [1 1 1;
    2 1 1;
    3 1 1;
    4 1 1;
    5 1 1];

lambda = 690;

nirsProbe = nirs.Probe( srcPos, detPos, link, lambda, [], [], srcDir, detDir );

% define optical properties
optProp(1) = nirs.OpticalProperties( .01, 1, 1.3, 690 );

% % % create forward model object
slabModel = nirs.SlabForwardModel( nirsProbe, optProp );
[slabJ, slabMeas] = slabModel.jacobian();

mcxModel = nirs.MCXForwardModel( nirsImage, optProp, nirsProbe );
mcxModel.nPhotons = 1e7;
mcxModel.nRepetitions = 1;
% mcxMeas = mcxModel.measurement();
[mcxJ, mcxMeas] = mcxModel.jacobian();

% measurement
% meas = fwdModel.measurement();
% [~,meas] = mcxModel.sensitivityMatrix();

% test with iterativeSolver
mcxModel.nPhotons = 1e7;
mcxModel.nRepetitions = 1;
mcxModel.optProp(1).mua = .02;
mcxModel.optProp(1).mus = 1.2;

optProp = nirs.mcxLayeredInverse( mcxMeas, mcxModel, @nirs.trustRegionKernel )


%%


clear

% segmented image
nirsImage = nirs.Image( uint8(ones(200,200,100)) );
nirsImage.volume(:,:,10:end) = 2;

% make a probe
detPos = [70 100 1];
srcPos = [90 100 1; 95 100 1; 100 100 1; 105 100 1; 110 100 1];
srcDir = repmat([0 0 1],[size(srcPos,1) 1]);
detDir = repmat([0 0 1],[size(detPos,1) 1]);

link = [1 1 1;
    2 1 1;
    3 1 1;
    4 1 1;
    5 1 1];

lambda = 690;

nirsProbe = nirs.Probe( srcPos, detPos, link, lambda, [], [], srcDir, detDir );

% define optical properties
optProp(1) = nirs.OpticalProperties( .01, 1, 1.3, 690 );
optProp(2) = nirs.OpticalProperties( .03, 1, 1.3, 690 );

mcxModel = nirs.MCXForwardModel( nirsImage, optProp, nirsProbe );
mcxModel.nPhotons = 1e7;
mcxModel.nRepetitions = 1;
% mcxMeas = mcxModel.measurement();
[~, mcxMeas] = mcxModel.jacobian();

% test with iterativeSolver
mcxModel.nPhotons = 1e7;
mcxModel.nRepetitions = 1;

slabSol = nirs.slabSolver( mcxMeas )

mcxModel.optProp(1).mua = slabSol.mua;
mcxModel.optProp(2).mua = slabSol.mua;

% mcxModel.optProp(1).mua = .018;
% mcxModel.optProp(1).mus = 1;
% 
% mcxModel.optProp(2).mua = .018;
% mcxModel.optProp(2).mus = 1;

optProp = nirs.mcxLayeredInverse( mcxMeas, mcxModel, @nirs.trustRegionKernel )