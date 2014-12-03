% clear
% 
% % segmented image
% img = nirs.Image( uint8(ones(150,150,100)) );
% 
% % make a probe
% detPos = [50 75 1];
% srcPos = [70 75 1; 75 75 1; 80 75 1; 85 75 1; 90 75 1];
% srcDir = repmat([0 0 1],[size(srcPos,1) 1]);
% detDir = repmat([0 0 1],[size(detPos,1) 1]);
% 
% link = [1 1 1;
%     2 1 1;
%     3 1 1;
%     4 1 1;
%     5 1 1;
%     1 1 2;
%     2 1 2;
%     3 1 2;
%     4 1 2;
%     5 1 2];
% 
% lambda = [690 830];
% 
% probe = nirs.Probe( srcPos, detPos, link, lambda, [], [], srcDir, detDir );
% 
% % optical properties
% tissue = nirs.OpticalTissue( 25e-6, 10e-6, 1, 1.3 );
% 
% % forward model
% fwdModel = nirs.LayeredSpectralModel( img, tissue, probe );
% fwdModel.fwdModel.nRepetitions = 10;
% 
% % measurement
% meas = fwdModel.measurement();
% 
% % fit slab
% nirs.spectralSlabSolver( meas )
% 
% % iterative fit
% invModel = fwdModel;
% invModel.optTissue = nirs.OpticalTissue( 35e-6, 15e-6, 1, 1.3 );
% invModel.fwdModel.nRepetitions = 5;
% nirs.mcxSpectralInverse( meas, invModel, @nirs.trustRegionKernel )

%%

clear

cd /home/barker
startup;

cd /home/barker/PhD_Data/+nirs/testScripts

% segmented image
R = 100;
vol = zeros((2*R+3)*[1 1 1]);

c = ceil(size(vol)/2);
[x,y,z] = meshgrid(...
    1:size(vol,1),...
    1:size(vol,2),...
    1:size(vol,3)...
    );

r = sqrt( (x(:)-c(1)).^2 ...
    + (y(:)-c(2)).^2 ...
    + (z(:)-c(3)).^2 ...
    );

vol( r < R ) = 1;
vol( r < (R-10) ) = 2;
vol( r < (R-12) ) = 3;

img = nirs.Image(uint8(vol));

% make a probe
srcPos = R*[cos(0/R)  sin(0/R) 0;
    cos(-2.5/R) sin(-2.5/R) 0];

detPos = R*[cos(20/R)  sin(20/R) 0;
    cos(25/R)  sin(25/R) 0;
    cos(30/R)  sin(30/R) 0;
    cos(35/R)  sin(35/R) 0];


srcPos(:,1) = srcPos(:,1) + c(1);
srcPos(:,2) = srcPos(:,2) + c(2);
srcPos(:,3) = srcPos(:,3) + c(3);

detPos(:,1) = detPos(:,1) + c(1);
detPos(:,2) = detPos(:,2) + c(2);
detPos(:,3) = detPos(:,3) + c(3);

srcDir = normr(repmat(c,[size(srcPos,1) 1]) - srcPos);
detDir = normr(repmat(c,[size(detPos,1) 1]) - detPos);

lambda = [690 830];

[X,Y,Z] = meshgrid(1:size(srcPos,1),1:size(detPos,1),1:length(lambda));

link = [X(:) Y(:) Z(:)];

probe = nirs.Probe( srcPos, detPos, link, lambda, [], [], srcDir, detDir );

% optical properties
tissue(1) = nirs.OpticalTissue( 60e-6, 40e-6, 1.2, 1.45 ); % skull 82/43 uM; skin 36/40 uM; mean 59/42
tissue(2) = nirs.OpticalTissue( 40e-6, 20e-6, 1.2, 1.45 ); % brain 70% /60 uM

% forward model
fwdModel = nirs.LayeredSpectralModel( img, tissue, probe );
fwdModel.fwdModel.nRepetitions = 1;

%% measurement
N = 16;

% meas = fwdModel.measurement();
% 
% for i = 2:N
%     meas(i) = fwdModel.measurement();
% end
% 
% save( 'test_data_8_distances_sphere.mat' )
 
load('test_data_8_distances_sphere.mat','meas')
for i = 2:N
    meas(1).data(i,:) = meas(i).data;
end

meas = meas(1);
meas.data = exp( mean( log(meas.data),1 ) );

%%
% iterative fit
invModel = fwdModel;
invModel.optTissue(1) = nirs.OpticalTissue( 50e-6, 50e-6, 1.2, 1.45 );
invModel.optTissue(2) = nirs.OpticalTissue( 40e-6, 40e-6, 1.2, 1.45 );
invModel.fwdModel.nRepetitions = 1;

fit = nirs.mcxSpectralInverse( meas, invModel, @nirs.fittingKernel);
% fit = nirs.mcxSpectralInverse( meas, invModel, @nirs.gradientDescentKernel );
% fit = nirs.mcxSpectralInverse( meas, invModel, @nirs.trustRegionKernel );
% fit = nirs.mcxSpectralInverse( meas, invModel, @nirs.gaussNewtonKernel_w_SVD );

save( 'testSpectral.mat' );