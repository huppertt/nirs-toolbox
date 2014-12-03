clear

cd /home/barker
startup;

cd /home/barker/PhD_Data/+nirs2/demo/

%%
img = nirs2.utilities.make_sphere( 100-cumsum([0 5 7 1.5]), 1.0 );

%%
R = 100;
d = (10:5:40)';
lambda = [690 830];

srcPos = R*[cos(0/R) sin(0/R) 0];
detPos = R*[cos(d/R) sin(d/R) zeros(size(d))];

srcDir = -normr(srcPos);
detDir = -normr(detPos);

[X,Y,Z] = meshgrid(1:size(srcPos,1),1:size(detPos,1),1:length(lambda));

link = [X(:) Y(:) Z(:)];

probe = nirs2.Probe( srcPos, detPos, link, lambda );
probe.srcDir = srcDir;
probe.detDir = detDir;

%%
so2 = 0.7; hbt = 60; lambda = [690 830];
prop{1} = nirs2.tissues.skin( lambda );
prop{2} = nirs2.tissues.bone( lambda );
prop{3} = nirs2.tissues.water( lambda ); % csf
prop{4} = nirs2.tissues.brain( so2, hbt, lambda );

%%
fwdModel = nirs2.MCXForwardModel( img, prop, probe, 110 );

[J,meas] = fwdModel.spectralJacobian();

fitProp = nirs2.solvers.slabSolver( meas );

save( 'sphere_demo.mat' )