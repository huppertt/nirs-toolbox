clear

cd /home/barker
startup;

cd /home/barker/PhD_Data/+nirs2/demo/

%%
% img = nirs2.utilities.make_sphere( 100-cumsum([0 5 5 1.5]), 1.0 );
img = nirs2.utilities.make_sphere( 100-cumsum([0 10]), 1.0 );

%%
R = 100;
s = 20*[0 -1/2]';
d = [10:20:70]';
lambda = [690 830];

srcPos = R*[cos(s/R) sin(s/R) zeros(size(s))];
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
prop{1} = nirs2.SpectralProperties( 0.5,100,lambda );
prop{2} = nirs2.tissues.brain( so2, hbt, lambda );

% prop{1} = nirs2.tissues.skin( lambda );
% prop{2} = nirs2.tissues.bone( lambda );
% prop{3} = nirs2.tissues.water( lambda ); % csf
% prop{4} = nirs2.tissues.brain( so2, hbt, lambda );

%%
fwdModel = nirs2.MCXForwardModel( img, prop, probe, 110 );

if ~exist('./2layer_sphere.mat','file')
    meas = fwdModel.measurement();
    for i = 2:16
        tmp = fwdModel.measurement();
        meas.data(i,:) = tmp.data;
    end
    save('./2layer_sphere.mat')
end

load('./2layer_sphere.mat','meas')

% meas.data = exp( mean(log(meas.data),1) );

%%
% img = nirs2.utilities.make_sphere( 100-cumsum([0 10 1.5]), 1.0 );
img = nirs2.utilities.make_sphere( 100-cumsum([0 10]), 1.0 );
 
% % superficial
% iProp{1} = nirs2.SpectralProperties( 0.6, 80, [690 830] );
% iProp{1}.mus = [2.3 1.7];
% % csf
% prop{2} = nirs2.tissues.water( lambda ); 
% % brain
% iProp{3} = nirs2.SpectralProperties( 0.7, 60, [690 830] );

iProp{1} = nirs2.SpectralProperties( 0.7, 100, [690 830] );
iProp{2} = nirs2.SpectralProperties( 0.7, 60, [690 830] );

% slabProp = nirs2.solvers.slabSolver( meas );
% iProp = {slabProp, slabProp};

invModel = nirs2.MCXForwardModel( img, iProp, probe, 110 );

fitProp = nirs2.solvers.spectralSolver( meas, invModel, logical([1 1]) );

save( 'sphere_fitting_demo.mat' )