%% Forward Model Demo 
%This script will demonstrate how to run the nirs.forward and nirs.inverse 
% solvers for the photon migration/image reconstruction problem
% This demo will show the intefaces to NIRFAST and NIRFASTBEM models for
% the forward model.

%% NIRFAST FEM model demo
% First, lets create a slab-model mesh
% nodes

xrange=[-80:5:80];
yrange=[-80:5:80];
zrange=[0:5:80];


[x, y, z] = ndgrid(xrange, yrange, zrange);
n = [x(:) y(:) z(:)];

% elements
e = delaunay(x(:), y(:), z(:));

% faces
f = [e(:, 1:3); e(:,[1 2 4]); e(:, [1 3 4]); e(:, 2:4)];
f = unique(sort(f,2), 'rows');

% This command will create a nirs.core.Mesh data type
mesh = nirs.core.Mesh(n, f, e);

%The mesh can be drawn using the command
mesh.draw();


% Next, let's create a probe for the simulation
srcPos = [-30 0 0];
detPos = [-20 0 0;
           -10 0 0;
           0 0 0;
           10 0 0;
           20 0 0;
           30 0 0];

link = [1 1 690;
    1 1 830;
    1 2 690;
    1 2 830;
    1 3 690;
    1 3 830;
    1 4 690;
    1 4 830;
    1 5 690;
    1 5 830;
    1 6 690;
    1 6 830;];

link = table(link(:,1), link(:,2), link(:,3), ...
    'VariableNames', {'source', 'detector', 'type'});
probe = nirs.core.Probe(srcPos, detPos, link);

% You can add "fiducial" points to a mesh which can be included while
% drawing.  The fiducials field needs to have a name, type, x,y,z, and a
% boolean Draw column (only true are drawn)

mesh.fiducials=[probe.optodes table(repmat(true,height(probe.optodes),1),'VariableNames',{'Draw'})];
mesh.transparency=.2;  % You can also adjust the transpency of the mesh for drawing
mesh.draw();

% The optical properties of the mesh are set using the nirs.media functions
% There are several predefined medias including brain, skull, csf  
% This defines mua/mus for brain at 70% O2 and 50 uM HbT
lambda=unique(probe.link.type)';  %List of wavelengths from the probe
prop{1} = nirs.media.tissues.brain(0.7, 50,lambda);

% each node is region 1
mesh.regions = ones(size(n,1), 1);

%If we had two regions, we would use
%>> prop{2} = nirs.media.tissues.skin(lambda);


% Now, create the forward model
fwdFEM = nirs.forward.NirfastFEM();
fwdFEM.mesh  = mesh;
fwdFEM.probe = probe;
fwdFEM.prop  = prop;
% The forward model class contains the methods for simulating measurements
% and computing the jacobian function.  Both CW and FD-NIRS are supports
% (but not time domain at this point).
% The forward model class exists for NIRFAST, NIRFAST-BEM, MCEXTREME, and the
% Semi-infinite slab models 

% To compute the measurement from this geometry
meas_FEM = fwdFEM.measurement();

% We can also compute the jacobian (in optical density)
[J, meas_FEM] = fwdFEM.jacobian();

% or the spectral jacobian (incorporates the beer-lambert law to give
% HbO2/Hb (and mus for the FD-NIRS models)
[J, meas_FEM] = fwdFEM.jacobian('spectral');


% Let's do the same thing using the NIRFAST-BEM model
%% NIRFAST BEM model demo
% First, lets create a slab-model mesh
% nodes
n=[];
[x,y,z]=ndgrid(xrange,yrange,zrange(1));
n = [n; x(:) y(:) z(:)];
[x,y,z]=ndgrid(xrange,yrange,zrange(end));
n = [n; x(:) y(:) z(:)];
[x,y,z]=ndgrid(xrange,yrange(1),zrange);
n = [n; x(:) y(:) z(:)];
[x,y,z]=ndgrid(xrange,yrange(end),zrange);
n = [n; x(:) y(:) z(:)];
[x,y,z]=ndgrid(xrange(1),yrange,zrange);
n = [n; x(:) y(:) z(:)];
[x,y,z]=ndgrid(xrange(end),yrange,zrange);
n = [n; x(:) y(:) z(:)];
n=unique(n,'rows');

% elements
e = delaunay(n);

% faces
f = [e(:, 1:3); e(:,[1 2 4]); e(:, [1 3 4]); e(:, 2:4)];
f = unique(sort(f,2), 'rows');

% This command will create a nirs.core.Mesh data type
mesh = nirs.core.Mesh(n, f);

fwdBEM = nirs.forward.NirfastBEM();
fwdBEM.mesh  = mesh;
% To add extra boundaries use:
% fwd.mesh(2)=mesh2, etc
fwdBEM.probe = probe;
fwdBEM.prop  = prop;

% To compute the measurement from this geometry
meas_BEM = fwdBEM.measurement();

% We can also compute the jacobian (in optical density)
%[J,meas_BEM] = fwdBEM.jacobian();

% or the spectral jacobian (incorporates the beer-lambert law to give
% HbO2/Hb (and mus for the FD-NIRS models)
%[J, meas_BEM] = fwdBEM.jacobian('spectral');

%% Now, let's try the slab model
fwdSlab = nirs.forward.SlabModel;
fwdSlab.Fm=0;
fwdSlab.probe=probe;
fwdSlab.prop=prop{1};

meas_Slab = fwdSlab.measurement();
%[J,meas_Slab]=fwdSlab.jacobian;

% Let's compare the models
figure; hold on;
plot(probe.distances,log(abs(meas_BEM.data)),'b');
plot(probe.distances,log(abs(meas_FEM.data)),'r');
plot(probe.distances,log(abs(meas_Slab.data)),'g');