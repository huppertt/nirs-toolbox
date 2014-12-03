clear

cd /home/barker
startup;

cd /home/barker/PhD_Data/+nirs2/demo/

%%
img = nirs2.Image( uint8(zeros(200,200,150)),0.5*[1 1 1],[100 100 5] );
img.vol(:,:,5:end) = 1;

%%
srcPos = [0 0 0];
detPos = [-35 0 0; 
    -25 0 0;
    -15 0 0;
    10 0 0;
    20 0 0;
    30 0 0];

% srcPos = srcPos + repmat( [100 100 5], [size(srcPos,1) 1] );
% detPos = detPos + repmat( [100 100 5], [size(detPos,1) 1] );

srcDir = repmat( [0 0 1],[size(srcPos,1) 1] );
detDir = repmat( [0 0 1],[size(detPos,1) 1] );

lambda = [690 830];

[X,Y,Z] = meshgrid( 1:size(srcPos,1), 1:size(detPos,1), 1:length(lambda) );
link = [X(:) Y(:) Z(:)];

probe = nirs2.Probe( srcPos, detPos, link, lambda );
probe.srcDir = srcDir;
probe.detDir = detDir;

%%
so2 = 0.7; hbt = 60;
prop = nirs2.tissues.brain( so2, hbt, lambda );

%%
Fm = 110; % modulation frequency
fwdModel = nirs2.MCXForwardModel( img, prop, probe, Fm );

meas = fwdModel.measurement();

%%

fitProp = nirs2.solvers.slabSolver( meas );

save( 'slab_demo.mat' )