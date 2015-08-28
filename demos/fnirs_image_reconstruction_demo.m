% Image reconstruction demo

%% First, lets set up the slab model for generating the forward model and the simulated data
xrange=[-100:5:100];
yrange=[-20:5:50];
zrange=[0:5:20];

[x, y, z] = ndgrid(xrange, yrange, zrange);
n = [x(:) y(:) z(:)];

% elements
e = delaunay(x(:), y(:), z(:));

% faces
f = [e(:, 1:3); e(:,[1 2 4]); e(:, [1 3 4]); e(:, 2:4)];
f = unique(sort(f,2), 'rows');

% This command will create a nirs.core.Mesh data type
mesh = nirs.core.Mesh(n, f, e);

%% Now let's create the probe
noise = nirs.testing.simARNoise();
stim  = nirs.testing.randStimDesign(noise.time, 10, 10, 2);

probe = noise.probe;
% This probe is already [-80 to 80mm, 0 to 25mm, 0] 

% If we add the probe as fiducial markers to the mesh, we can view it
mesh.fiducials=[probe.optodes table(repmat(true,height(probe.optodes),1),'VariableNames',{'Draw'})];
mesh.transparency=.2;  % You can also adjust the transpency of the mesh for drawing
mesh.draw();


% Let's put a target image in the mesh

targetLoc=[20 0 10; -20 5 10];
extent=[15 15];
Amp=[10 -2; 10 -3];  % Amplitude of the two conditions x HbO2/Hb


mesh.regions=ones(size(mesh.nodes,1),1);
lambda=unique(probe.link.type)';  %List of wavelengths from the probe
prop{1} = nirs.media.tissues.brain(0.7, 50,lambda);

% Now, create the forward model
fwdFEM = nirs.forward.NirfastFEM();
fwdFEM.mesh  = mesh;
fwdFEM.probe = probe;
fwdFEM.prop  = prop;
[J, meas_FEM_baseline] = fwdFEM.jacobian('spectral');


for idx=1:size(targetLoc,1)
    mask(:,idx)=[Amp(idx,1)*(sqrt(sum((mesh.nodes-...
        ones(size(mesh.nodes,1),1)*targetLoc(idx,:)).^2,2))<=extent(idx)); ...
        Amp(idx,2)*(sqrt(sum((mesh.nodes-...
        ones(size(mesh.nodes,1),1)*targetLoc(idx,:)).^2,2))<=extent(idx))];
end

% to draw use:
%>> mesh.draw(mask(1:end/2,1))

[data,truth] = nirs.testing.simDataImage(fwdFEM, noise, stim, mask);

[data,truth] = nirs.testing.simDataImage(fwdFEM, noise, stim);

%% Now solve the model
j = nirs.modules.TrimBaseline();
j = nirs.modules.OpticalDensity(j);
j = nirs.modules.Resample(j);
j = nirs.modules.AR_IRLS(j);

SubjStats = j.run(data);

j = nirs.modules.ImageReconMFX();

% The jacobian is a dictionary indexed by the subject name (or default)
j.jacobian('default')=J;
j.probe('default')=noise.probe;
j.mesh=mesh;  
j.formula = 'beta ~ -1 + cond';  % Simple fixed effects model

j.basis=nirs.inverse.basis.identity(mesh);

%j.basis=nirs.inverse.basis.gaussian(mesh,.1);
% Now create the priors in the model
prior.hbo(:,1)=zeros(size(J.hbo,2),1);
prior.hbr(:,1)=zeros(size(J.hbr,2),1);
% 
%  prior.hbo(:,2)=mask(1:end/2,1);
%  prior.hbr(:,2)=mask(end/2+1:end,1);
% % 
% % prior.hbo(:,3)=mask(1:end/2,2);
% % prior.hbr(:,3)=mask(end/2+1:end,2);
 prior2.hbo(:,1)=zeros(size(J.hbo,2),1);
 prior2.hbr(:,1)=zeros(size(J.hbr,2),1);
% % 
%  prior2.hbo(:,2)=mask(1:end/2,2);
%  prior2.hbr(:,2)=mask(end/2+1:end,2);
% The fields and dimensions in in the prior need to match that of the Jacobian
j.prior('A')=prior;
j.prior('B')=prior2;
% Prior is a dictionary and uses the names of the stimulus conditions in the model

ImageStats=j.run(SubjStats);

ImageStats.draw('tstat',[-5 5])


ROCtest=nirs.testing.ChannelStatsROC;
ROCtest.simfunc=@()nirs.testing.simDataImage(fwdFEM, noise, stim);
j = nirs.modules.TrimBaseline();
j = nirs.modules.OpticalDensity(j);
j = nirs.modules.Resample(j);
j = nirs.modules.AR_IRLS(j);
j = nirs.modules.ImageReconMFX(j);
j.jacobian('default')=J;
j.probe('default')=noise.probe;
j.mesh=mesh;  
j.formula = 'beta ~ -1 + cond';  % Simple fixed effects model
j.basis=nirs.inverse.basis.identity(mesh);
prior.hbo(:,1)=zeros(size(J.hbo,2),1);
prior.hbr(:,1)=zeros(size(J.hbr,2),1);
j.prior('default')=prior;

num_iter = 10;  

% Assign the job to the ROCtest pipeline
ROCtest.pipeline=j;

% FInally run the ROC for # of iterations
ROCtest=ROCtest.run(num_iter);

ROCtest.draw;
