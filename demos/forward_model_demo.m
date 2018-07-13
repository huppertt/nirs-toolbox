%% Forward Model Demo 
%This script will demonstrate how to run the nirs.forward and nirs.inverse 
% solvers for the photon migration/image reconstruction problem
% This demo will show the intefaces to NIRFAST and NIRFASTBEM models for
% the forward model.

%% NIRFAST FEM model demo
% First, lets create a slab-model mesh

Slab = nirs.core.Image;
Slab.dim = [5 5 5];
Slab.origin = [-100 -100 0];
Slab.description = 'Slab model for FEM/BEM models';
Slab.vol = ones(41,41,11);  % SLab from -100:5:100, -100:5:100, 0:5:50

% This command will create a nirs.core.Mesh data type
% This requires iso2mesh to be installed.  
% Cite:
% Qianqian Fang and David Boas, "Tetrahedral mesh generation from volumetric binary
% and gray-scale images," Proceedings of IEEE International Symposium on Biomedical 
% Imaging 2009, pp. 1142-1145, 2009
%               
% http://iso2mesh.sourceforge.net

% If the code is not installed, you will get an error and instructions on
% downloading from the next line.  Iso2Mesh is not distributed with the
% nirs-toolbox and needs to be seperately downloaded
      
mesh = Slab.convertFEMmesh();   % This will convert the (binary) Image/Vol to a mesh 
mesh = mesh.reducemesh(.2);  % This will resample the mesh using Iso2Mesh


%The mesh can be drawn using the command
% figure;
% mesh.draw();


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
mesh.regions = ones(size(mesh.nodes,1), 1);

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
%[J, meas_FEM] = fwdFEM.jacobian();

% or the spectral jacobian (incorporates the beer-lambert law to give
% HbO2/Hb (and mus for the FD-NIRS models)
%[J, meas_FEM] = fwdFEM.jacobian('spectral');


% Let's do the same thing using the NIRFAST-BEM model
%% NIRFAST BEM model demo
mesh = Slab.convertBEMmesh();   % This will convert the (binary) Image/Vol to a surface mesh 
mesh = mesh.reducemesh(.2);  % This will resample the mesh using Iso2Mesh

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


%% Now the MCextreme model

if(0)
    %This will run the MCextreme model but will only work if the comupter
    %supports GPU processing and has the appropriate Nvidia CUDA drivers
    %installed
    fwdMCX = nirs.forward.MCXLab;
    fwdMCX.probe=probe;
    fwdMCX.prop=prop;
    
    fwdMCX.image=Slab;
    fwdMCX.Fm=0;
    
    meas_Slab = fwdMCX.measurement();
    
end

% The theoretical model from
% Fantini et al 1994 Applied Optics 33(22) 5204-5213
r=probe.distances;

[~,iLambda]=ismember(probe.link.type,prop{1}.lambda);
mua=prop{1}.mua(iLambda)';
kappa=prop{1}.kappa(iLambda)';
Udc = 1./r.*exp(-r.*sqrt(mua./kappa));


% Let's compare the models

figure; hold on;
lst1=find(iLambda==1);  % First wavelength
plot(probe.distances(lst1),log10(abs(meas_BEM.data(lst1))),'b');
plot(probe.distances(lst1),log10(abs(meas_FEM.data(lst1))),'r');
plot(probe.distances(lst1),log10(abs(meas_Slab.data(lst1))),'g');
plot(probe.distances(lst1),log10(Udc(lst1)),'k');

lst2=find(iLambda==2);  % Second wavelength
plot(probe.distances(lst2),log10(abs(meas_BEM.data(lst2))),'b--');
plot(probe.distances(lst2),log10(abs(meas_FEM.data(lst2))),'r--');
plot(probe.distances(lst2),log10(abs(meas_Slab.data(lst2))),'g--');
plot(probe.distances(lst2),log10(Udc(lst2)),'k--');

legend({'BEM Model','FEM Model','Slab Model','Analytic'})
xlabel('Src-Det distance (mm)');
ylabel('Log10(Y)');
