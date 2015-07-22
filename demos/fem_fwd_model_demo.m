%% mesh
% nodes
[x, y, z] = ndgrid(-50:10:50, -50:10:50, 0:10:50);

n = [x(:) y(:) z(:)];

% elements
e = delaunay(x(:), y(:), z(:));

% faces
f = [e(:, 1:3); e(:,[1 2 4]); e(:, [1 3 4]); e(:, 2:4)];
f = unique(sort(f,2), 'rows');

% mesh
mesh = nirs.core.Mesh(n, f, e);

%% probe
srcPos = [0 0 0];
detPos = [-20 0 0;
        20 0 0];

link = [1 1 690;
    1 1 830;
    1 2 690;
    1 2 830];

link = table(link(:,1), link(:,2), link(:,3), ...
    'VariableNames', {'source', 'detector', 'type'});

probe = nirs.core.Probe(srcPos, detPos, link);

%% properties
% brain at 70% O2 and 50 uM HbT
prop{1} = nirs.media.tissues.brain(0.7, 50);

% each node is region 1
mesh.regions = ones(size(n,1), 1);

%% forward model
fwd = nirs.forward.NirfastFEM();
fwd.mesh  = mesh;
fwd.probe = probe;
fwd.prop  = prop;

% measurements
meas = fwd.measurement();

% standard jacobian
[J, meas] = fwd.jacobian();

% spectral jacobian
[J, meas] = fwd.jacobian('spectral');

